# Handle access to optking module

# Some handy constants
BEND = "bend"
STRE = "stre"
LIN1 = "lin1"
LIN2 = "lin2"
TORS = "tors"
OUT  = "out"

module Psi
  class OptKing
        
    # Mixin the InputGenerator
    include InputGenerator
    # We don't mix in the Executor because optking returns error codes corresponding to certain actions
    
    def initialize(task_obj)
      set_task task_obj
    end
    
    # Generate the :fixed_info tag for the input file
    def generate_fixedintco(input_file, fixedintco)
      if fixedintco != nil
        input_file.puts "fixed_intco:("
        
        # Go through fixedintco and generate the input
        sub_result = generate_value(fixedintco)
        input_file.puts(sub_result.string)
        
        input_file.puts ")"
      end
      return input_file
    end
    
    def execute(input_file, binary)
      puts "DEBUG: For #{binary} would have used:\n#{input_file.string}" if Psi::check_commands == true
      prefix = ""
      prefix = "-p #{get_task.prefix}" if get_task.prefix != nil
      
      if Psi::check_commands != true
        # Run optking with our input file
        if Psi::quiet == false
          `echo "#{input_file.string}" | #{binary} #{prefix} -f -`
        else
          `echo "#{input_file.string}" | #{binary} #{prefix} -f - >& /dev/null`
        end
        if $?.exited?
          # Return the status code to the caller so that it can handle what to do next
          return $?.exitstatus
        else
          puts "Error occurred in OptKing::execute."
          exit 1
        end
      end
    end
  end
  
  # Add optking ability to Task
  class Task
    
    def optking(*args)
      
      # The user MUST provide a code block
      if block_given? == false
        puts red("Error: ") + "You must provide a code block for optking."
        exit 1
      end
      
      # Convert to a hash
      args_hash = args[0]
      # Initial values
      result = 0
      max_geometry_cycles = 20
      current_cycle = 0
      save_energy = false
      fixedintco = nil
      cleanup = true

      # By default, we do by energy points unless the requested wavefunction supports analytic gradients
      set_gradients(false)
      
      # Check to see what we support:
      if wavefunction.casecmp("scf") == 0 and SCF.supports_analytic_gradients[reference] == true
        set_gradients(true)
      elsif wavefunction.casecmp("ccsd") == 0 and CCEnergy.supports_analytic_gradients[reference] == true
        set_gradients(true)
      elsif wavefunction.casecmp("ccsd_t") == 0 and CCTriples.supports_analytic_gradients[reference] == true
        set_gradients(true)
      elsif wavefunction.casecmp("mp2") == 0 and MP2.supports_analytic_gradients[reference] == true
        set_gradients(true)
      elsif wavefunction.casecmp("detci") == 0 and DetCI.supports_analytic_gradients[reference] == true
        set_gradients(true)
      end
      
      # Check to see if the user requested gradients, if they did and get_gradients is false
      # warn the user and do by energy points
      if args_hash != nil
        # Did the user specify fixed coordinates?
        if args_hash.has_key?(:fixed_intco)
          puts "Doing a constrained geometry optimization."
          fixedintco = args_hash[:fixed_intco]
        end
        
        if args_hash.has_key?(:cleanup)
          cleanup = args_hash[:cleanup]
        end
        
        if args_hash.has_key?(:gradients) and args_hash[:gradients] == true and get_gradients == false
          puts red("WARNING: ") + "Analytic gradients requested, but are not supported."
          puts "         Switching to gradients by energy points."
        end

        # Did the user explicitly request by energy points?
        if args_hash.has_key?(:gradients) and args_hash[:gradients] == false
          # User requested NO gradients.
          set_gradients(false)
        end
        
        # Does user want to save result of code block to checkpoint for optking to use?
        if args_hash.has_key?(:save_energy)
          save_energy = args_hash[:save_energy]
        end
      end
      
      # Did user request a limited number of optimization cycles?
      if args_hash != nil
        if args_hash.has_key?(:max_cycles)
          max_geometry_cycles = args_hash[:max_cycles]
        end
      end
      
      # We are done with :gradients option, drop it from the list so that it isn't passed to optking
      args_hash.delete(:gradients) unless args_hash == nil
      args_hash.delete(:max_cycles) unless args_hash == nil
      args_hash.delete(:fixed_intco) unless args_hash == nil
      
      # Report what we are doing
      if get_gradients == true
        puts "Optimizing geometry via analytic gradients.\n"
      else
        puts "Optimizing geometry via energy points.\n"
      end
      
      # Create a new optking object
      optking_obj = Psi::OptKing.new self

      # Form the input hash and generate input file
      input_hash = { "reference" => reference, "wfn" => wavefunction }
      
      # Make sure gradients are turned off or on
      if get_gradients == false
        input_hash["dertype"] = "none"
      else
        input_hash["dertype"] = "first"
      end
      
      # We are doing optimization
      input_hash["jobtype"] = "opt"
      input_hash = input_hash.merge(args_hash) unless args_hash == nil

      # Handle the geometry
      if @zmat != nil and @geometry == nil
        input_hash["zmat"] = @zmat
      elsif @zmat == nil and geometry != nil
        input_hash["geometry"] = @geometry
      elsif @zmat == nil and @geometry == nil
        puts red("Error: ") + "Neither geometry nor zmat are set."
        exit 1
      else
        puts red("Error: ") + "Both geometry and zmat are set. One must be nil."
        exit 1
      end
      
      # Generate optking input file
      input_file = optking_obj.generate_input input_hash
      
      # If using fixedintco need to add it to the input file
      input_file = optking_obj.generate_fixedintco(input_file, fixedintco) unless fixedintco == nil
            
      puts "Optimizing geometry..."
      Psi::indent_puts

      # REQUIRED: Run input to generate initial checkpoint file contents
      input(args_hash)

      # OPTIONAL: Run a single point so that optking has something to work with if we are doing energy points
      if get_gradients == false        
        loop do
          puts "Computing reference geometry and energy:"
          Psi::indent_puts
          result = yield(input_hash)

          # Save energy to checkpoint if needed
          self.etot = result if save_energy == true
          Psi::unindent_puts

          # Single point obtained, form displacements
          current_cycle += 1
          number_displacements = optking_obj.execute(input_file, Psi::Commands::FINDIF_DISP_SYMM)
          if number_displacements == 1
            # Something must be wrong
            puts "Error: Unable to generate displacements."
            exit 1
          end
          puts "Cycle #{current_cycle}: computing energies for #{number_displacements} displacements:"
          Psi::indent_puts
          
          # Attempt to do the displacements
          number_displacements.times do |i|
            # Do pre-displacement commands
            optking_obj.execute(input_file, Psi::Commands::FINDIF_NEXT)
            input(:binary => Psi::Commands::FINDIF_INPUT, :quiet => true)

            puts "Displacement #{i+1}"
            Psi::indent_puts

            # Get the energy from user
            result = yield(input_hash)

            # Save energy to checkpoint if needed
            self.etot = result if save_energy == true
            Psi::unindent_puts

            # Do post-displacement commands
            optking_obj.execute(input_file, Psi::Commands::FINDIF_ENERGY_SAVE)
          end
          Psi::unindent_puts
          
          # Finished with displacements
          # Tell optking to compute gradients from energies
          result = optking_obj.execute(input_file, Psi::Commands::FINDIF_GRAD_ENERGY)
          # Tell optking to take a step.
          # If result == 2 geometry is optimized
          # If result == 0 geometry is not optimized
          # If result == 1 optking failure
          result = optking_obj.execute(input_file, Psi::Commands::GEOMUPDATE)
          break if result == 2 or result == 1
          
          # Is it time to stop?
          if current_cycle >= max_geometry_cycles
            result = 3
            break
          end
        end
        
        # Report results
        if result == 2
          puts "Optimization complete."
        elsif result == 3
          puts "Max cycles reached. Stopping optimization."
        else
          puts "Error during optimization."
          exit 1
        end
      else
        # Optimizations by analytic gradients are much easier
        loop do
          current_cycle += 1
          puts "Cycle #{current_cycle}:"
          Psi::indent_puts
          
          # Yield to the user code block passing the need input_hash
          result = yield(input_hash)
          
          # Update geometry
          result = optking_obj.execute(input_file, Psi::Commands::GEOMUPDATE)
          break if result == 2 or result == 1
          
          # Is it time to stop?
          if current_cycle >= max_geometry_cycles
            result = 3
            break
          end
          
          Psi::unindent_puts
        end
        
        # Report results
        if result == 2
          puts "Optimization complete."
        elsif result == 3
          puts "Max cycles reached. Stopping optimization."
        else
          puts "Error during optimization."
          exit 1
        end
      end
      
      Psi::unindent_puts
      Psi::unindent_puts # I'm missing one somewhere
      
      # Tell methods that follow that an optimization was completed.
      set_optimization_complete true
      
      # Unless told not to, run psiclean
      clean if cleanup == true
      
      # If the user gave a zmatrix read in the new zmatrix from checkpoint
      if @zmat != nil
        chkpt_zmat = ZEntry.new self
        zmat_arr = chkpt_zmat.to_a
        
        # zmat_arr should "look" like @zmat, go through it and update the variables
        # if a variable is a symbol then update the global variable referred to by the symbol
        count = zmat_arr.length - 1
        0.upto(count) do |x|
          if @zmat[x].length >= 3
            if @zmat[x][2].kind_of?(Symbol)
              # Set the global variable that is named from the symbol to the zmat_arr value
              eval("$" + @zmat[x][2].id2name + "=zmat_arr[x][2]")
            else
              @zmat[x][2] = zmat_arr[x][2]
            end
          end
          if @zmat[x].length >= 5
            if @zmat[x][4].kind_of?(Symbol)
              # Set the global variable that is named from the symbol to the zmat_arr value
              eval("$" + @zmat[x][4].id2name + "=zmat_arr[x][4]")
            else
              @zmat[x][4] = zmat_arr[x][4]
            end
          end
          if @zmat[x].length >= 7
            if @zmat[x][6].kind_of?(Symbol)
              # Set the global variable that is named from the symbol to the zmat_arr value
              eval("$" + @zmat[x][6].id2name + "=zmat_arr[x][6]")
            else
              @zmat[x][6] = zmat_arr[x][6]
            end
          end
        end
      end
      
      # Return the new total energy of the system
      etot
    end
  end
end

# User is expected to provide a code block that tells optking what to do
# for a single step.
# Accepts the following options:
#  :gradients => true/false      default depends on wavefunction, user can only override with false
#  :max_cycles => number         default 20
#  :save_energy => true/false    default false (optking reads energy from checkpoint), true (save result of
#                                  code block to checkpoint which optking will read)
def optking(*args)
  # Convert to a hash
  args_hash = args[0]
  Psi::global_task.optking(args_hash) do |input_hash|
    yield(input_hash)
  end
end
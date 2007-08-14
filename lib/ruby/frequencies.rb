module Psi
  class Frequencies
    # Mixin the InputGenerator
    include InputGenerator
    
    def initialize(task_obj)
      set_task task_obj
      # Override the psi_module_name
      set_psi_module_name "OptKing"
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
          puts "Error occurred in Frequencies::execute."
          exit 1
        end
      end
    end
  end
  
  # Add frequencies ability to Task
  class Task
    def frequencies(*args)
      # The user MUST provide a code block
      if block_given? == false
        puts red("Error: ") + "You must provide a code block for frequencies."
        exit 1
      end
      
      # Convert to a hash
      args_hash = args[0]
      
      # By default do frequencies by energy points unless the requested wavefunction
      # supports either second derivatives or first derivatives
      set_gradients(false)
      set_second_deriv(false)
      
      # Check to see what we support
      # First do gradients
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
      
      # Next check for second derivatives
      # Check to see what we support:
      if wavefunction.casecmp("scf") == 0 and SCF.supports_analytic_second_derivatives[reference] == true
        set_second_deriv(true)
      elsif wavefunction.casecmp("ccsd") == 0 and CCEnergy.supports_analytic_second_derivatives[reference] == true
        set_second_deriv(true)
      elsif wavefunction.casecmp("ccsd_t") == 0 and CCTriples.supports_analytic_second_derivatives[reference] == true
        set_second_deriv(true)
      elsif wavefunction.casecmp("mp2") == 0 and MP2.supports_analytic_second_derivatives[reference] == true
        set_second_deriv(true)
      elsif wavefunction.casecmp("detci") == 0 and DetCI.supports_analytic_second_derivatives[reference] == true
        set_second_deriv(true)
      end
      
      # Check for user arguments
      if args_hash != nil
        # Did the user request gradients?
        if args_hash.has_key?(:gradients) and args_hash[:gradients] == true and get_gradients == false
          puts red("WARNING: ") + "Analytic gradients requested, but are not supported."
          puts "         Switching to frequencies by energy points."
        elsif args_hash.has_key?(:gradients) and args_hash[:gradients] == true and get_second_deriv == true
          puts red("WARNING: ") + "Analytic second derivatives available, but analytic gradients requested."
          puts "         Switching to frequencies by analytic gradients."
          # If we have 2nd derivatives then we must have 1st derivatives
          set_second_deriv(false)
        elsif args_hash.has_key?(:gradients) and args_hash[:gradients] == false
          set_gradients(false)
          set_second_deriv(false)
        end
      end
      args_hash.delete(:gradients) unless args_hash == nil
      
      # Report what we are doing
      if get_second_deriv == true
        set_gradients(false)
        puts "\nComputing frequencies via analytic second derivatives."
      elsif get_gradients == true
        puts "\nComputing frequencies via finite differences of gradients."
      else
        puts "\nComputing frequencies via finite differences of energies."
      end
      
      # Create a new frequencies object
      frequencies_obj = Psi::Frequencies.new self
      
      # Form the input hash and generate input file.
      # Form the input hash and generate input file
      input_hash = { "reference" => reference, "wfn" => wavefunction }
      
      # Make sure derivatives are turned off or on
      if get_second_deriv == true
        input_hash["dertype"] = "second"
      elsif get_gradients == true
        input_hash["dertype"] = "first"
      else
        input_hash["dertype"] = "none"
      end
      
      # Handle the geometry
      if get_optimization_complete == false
        if @zmat != nil and @geometry == nil
          input_hash["zmat"] = @zmat
        elsif @zmat == nil and @geometry != nil
          input_hash["geometry"] = @geometry
        elsif @zmat == nil and @geometry == nil
          puts red("Error: ") +"Neither geometry nor zmat are set."
          exit 1
        else
          puts red("Error: ") + "Both geometry and zmat are set. One must be nil."
          exit 1
        end
      end
      
      # We are doing frequiencies
      input_hash["jobtype"] = "freq"
      input_hash = input_hash.merge(args_hash) unless args_hash == nil
      
      # Generate the input file
      input_file = frequencies_obj.generate_input input_hash
      
      puts "Computing frequencies..."
      Psi::indent_puts
      
      # Run input to generate initial checkpoint file ONLY if get_optimization_completed == false
      if get_optimization_complete == false
        input(args_hash)
      else
        # else run input with --chkptgeom
        puts "Using optimized geometry from checkpoint."
        #input(:binary => Psi::Commands::INPUTCHKPT, :quiet => true)
      end
      
      if get_second_deriv == true
        # Simply yield out to the code block, it should know what to do.
        yield(nil)
      elsif get_gradients == true
        # If geometry optimization was done then we do not need to do a single point energy
        if get_optimization_complete == false
          # Temporarily turn off gradients
          set_gradients(false)
          # Do the single point energy
          yield(nil)
          # Turn gradients back on
          set_gradients(true)
        end
        
        # Figure out the displacements needed for frequencies by gradients
        number_displacements = frequencies_obj.execute(input_file, Psi::Commands::FINDIF_DISP_FREQ_GRAD_CART)
        puts "Computing #{number_displacements} displacements."
        Psi::indent_puts
        
        # Attempt to do the displacements
        number_displacements.times do |i|
          # Do pre-displacement commands
          frequencies_obj.execute(input_file, Psi::Commands::FINDIF_NEXT)
          input(:binary => Psi::Commands::FINDIF_INPUT, :quiet => true)
          
          puts "Displacement #{i+1}"
          # Compute the gradient at this point
          Psi::indent_puts
          yield(nil)
          Psi::unindent_puts
          
          # Do post-displacement commands
          frequencies_obj.execute(input_file, Psi::Commands::FINDIF_GRAD_SAVE)
        end
        Psi::unindent_puts
        
        # Finished with displacements
        # Tell optking to compute frequencies from the gradients
        frequencies_obj.execute(input_file, Psi::Commands::FINDIF_FREQ_GRAD_CART)
      else
        # If geometry optimization was done then we do not need to do a single point energy
        if get_optimization_complete == false
          # Do the single point energy
          yield(nil)
        end
        
        # Figure out the displacements needed for frequencies by energy
        number_displacements = frequencies_obj.execute(input_file, Psi::Commands::FINDIF_DISP_FREQ_ENERGY_CART)
        puts "Computing #{number_displacements} displacements."
        Psi::indent_puts
        
        # Attempt to do the displacements
        number_displacements.times do |i|
          # Do pre-displacement commands
          frequencies_obj.execute(input_file, Psi::Commands::FINDIF_NEXT)
          input(:binary => Psi::Commands::FINDIF_INPUT, :quiet => true)
          
          puts "Displacement #{i+1}"
          # Compute the energy at this point
          Psi::indent_puts
          yield(nil)
          Psi::unindent_puts
          
          # Do post-displacement commands
          frequencies_obj.execute(input_file, Psi::Commands::FINDIF_ENERGY_SAVE)
        end
        Psi::unindent_puts
        
        # Finished with displacements
        # Tell optking to compute frequencies from the energies
        frequencies_obj.execute(input_file, Psi::Commands::FINDIF_FREQ_ENERGY_CART)
      end
      Psi::unindent_puts
    end
  end
end

def frequencies(*args)
  # Convert to a hash
  args_hash = args[0]
  Psi::global_task.frequencies(args_hash) do |input_hash|
    yield(input_hash)
  end
end
alias freq frequencies

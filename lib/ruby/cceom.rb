# Handle access to Psi ccenergy module

module Psi
  class CCEom
    
    # Mixin the InputGenerator
    include InputGenerator
    include Executor
    # Which references does this method support
    def self.supports_analytic_gradients
      { "rhf" => false, "rohf" => false, "uhf" => false, "twocon" => false }
    end
    def self.supports_analytic_second_derivatives
      { "rhf" => false, "rohf" => false, "uhf" => false, "twocon" => false }
    end
    
    def initialize(task_obj)
      @task = task_obj
      # Set the generic command for this class
      set_binary_command Psi::Commands::CCEOM
    end
  end
  
  # Add ccenergy ability to Task
  class Task
    def eom(*args)
      # convert to a hash
      args_hash = args[0]

      # Call transqt and ccsort unless told not to
      if args_hash != nil
        if args_hash.has_key?(:cchbar) == false or 
           (args_hash.has_key?(:cchbar) and args_hash[:cchbar] == true)
          cchbar(args_hash)
        end
      else
        cchbar(args_hash)
      end
      args_hash.delete(:cchbar) unless args_hash == nil

      # Create a new scf object
      eom_obj = Psi::CCEom.new self

      # Form the input hash and generate the input file
      input_hash = { }

      # Check to see if the function arguments have the reference, if so use it, otherwise use
      # global setting
      if args_hash == nil or args_hash.has_key?("reference") == false
        input_hash["reference"] = reference
      end

      # If we are doing analytic gradients make sure ccenergy knows
      if get_gradients == true and CCEom.supports_analytic_gradients[reference] == true
        input_hash["dertype"] = "first"
      else
        input_hash["dertype"] = "none"
      end
      
      # Check the wavefunction
      if args_hash == nil or args_hash.has_key?("wfn") == false
        input_hash["wfn"] = wavefunction
      end

      # Merge what we've done with what the user wants, user settings should override
      input_hash = input_hash.merge(args_hash) unless args_hash == nil

      # Run the ccenergy module, sending the input file as keyboard input
      puts "eom"
      eom_obj.execute(input_hash)
      
      # Check to see if we are supposed to compute analytical gradients, only do if
      #  1. It is supported
      #  2. The user requested it
      #  3. The wavefunction is CCSD
      if get_gradients == true and CCEom.supports_analytic_gradients[reference] == true and wavefunction.casecmp("eom_ccsd") == 0
        puts "gradient:"
        Psi::indent_puts
        cclambda(input_hash, :binary => Psi::Commands::CCLAMBDA_EXCITED)
        ccdensity(input_hash, :binary => Psi::Commands::CCDENSITY_CALC_XI)
        cclambda(input_hash, :binary => Psi::Commands::CCLAMBDA_ZETA)
        ccdensity(input_hash, :binary => Psi::Commands::CCDENSITY_USE_ZETA)
        oeprop(input_hash)
        
        # Tell transqt to do a back transformation
        transqt_input_hash = input_hash
        transqt_input_hash[:backtr] = true
        transqt(transqt_input_hash)
        deriv(input_hash)
        Psi::unindent_puts
      end
      
      eom_states_energy
    end
  end
end

# Create some global functions
# User can send additional input parameters to the function
def eom(*args)
  # convert to a hash
  args_hash = args[0]
  Psi::global_task.eom(args_hash)
end

# Handle access to Psi detci module

module Psi
  class DetCI
    
    # Mixin the InputGenerator
    include InputGenerator
    include Executor
    # Which references does this method support
    def self.supports_analytic_gradients
      { "rhf" => false, "rohf" => false, "uhf" => false, "twocon" => false }
    end
    
    def initialize(task_obj)
      @task = task_obj
      # Set the generic command for this class
      set_binary_command Psi::Commands::DETCI
    end
  end
  
  # Add detci ability to Task
  class Task
    def detci(*args)
      # convert to a hash
      args_hash = args[0]

      # Call transqt unless told not to
      if args_hash != nil
        if args_hash.has_key?(:transqt) == false or 
           (args_hash.has_key?(:transqt) and args_hash[:transqt] == true)
          transqt(args_hash)
        end
      else
        transqt(args_hash)
      end
      args_hash.delete(:transqt) unless args_hash == nil

      # Create a new scf object
      detci_obj = Psi::DetCI.new self

      # Form the input hash and generate the input file
      input_hash = { }

      # Check to see if the function arguments have the reference, if so use it, otherwise use
      # global setting
      if args_hash == nil or args_hash.has_key?("reference") == false
        input_hash["reference"] = reference
      end

      # If we are doing analytic gradients make sure detci knows
      if get_gradients == true and DetCI.supports_analytic_gradients[reference] == true
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
      puts "detci"
      detci_obj.execute(input_hash)
      
      etot
    end
  end
end

# Create some global functions
# User can send additional input parameters to the function
def detci(*args)
  # convert to a hash
  args_hash = args[0]
  Psi::global_task.detci(args_hash)
end

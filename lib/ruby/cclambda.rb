# Handle access to Psi cclambda module

module Psi
  class CCLambda
    # Mixin the InputGenerator
    include InputGenerator
    include Executor
    
    def initialize(task_obj, binary=nil)
      @task = task_obj
      # Set the generic command for this class
      if binary == nil
        set_binary_command Psi::Commands::CCLAMBDA
      else
        set_binary_command binary
      end
    end
  end
  
  # Add cclambda ability to Task
  class Task
    def cclambda(*args)
      # convert to a hash
      args_hash = args[0]

      binary = nil
      if args_hash != nil
        if args_hash.has_key?(:binary)
          binary = args_hash[:binary]
        end
      end
      
      # Create a new cclambda object
      cclambda_obj = Psi::CCLambda.new(self, binary)

      # Form the input hash and generate the input file
      input_hash = { }

      # Check to see if the function arguments have the reference, if so use it, otherwise use
      # global setting
      if args_hash == nil or args_hash.has_key?("reference") == false
        input_hash["reference"] = reference
      end

      # If we are doing analytic gradients make sure cchbar knows
      if get_gradients == true
        input_hash["dertype"] = "first"
      else
        input_hash["dertype"] = "none"
      end
      
      # Check the wavefunction
      if args_hash == nil or args_hash.has_key?("wfn") == false
        input_hash["wfn"] = wavefunction
      end

      # Merge what we've done with what the user wants
      input_hash = input_hash.merge(args_hash) unless args_hash == nil

      # Run the ccenergy module, sending the input file as keyboard input
      puts "cclambda"
      cclambda_obj.execute(input_hash)
    end
  end
end

# Create some global functions
# User can send additional input parameters to the function
def cclambda(*args)
  # convert to a hash
  args_hash = args[0]
  Psi::global_task.cclambda(args_hash)
end

# Handle access to Psi oeprop module

module Psi
  class OEProp
    # Mixin the InputGenerator
    include InputGenerator
    include Executor
    
    def initialize(task_obj)
      @task = task_obj
      # Set the generic command for this class
      set_binary_command Psi::Commands::OEPROP
    end
  end
  
  # Add oeprop ability to Task
  class Task
    def oeprop(*args)
      # convert to a hash
      args_hash = args[0]

      # Create a new cclambda object
      oeprop_obj = Psi::OEProp.new self

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
      puts "oeprop"
      oeprop_obj.execute(input_hash)
    end
  end
end

# Create some global functions
# User can send additional input parameters to the function
def oeprop(*args)
  # convert to a hash
  args_hash = args[0]
  Psi::global_task.oeprop(args_hash)
end

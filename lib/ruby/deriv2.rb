
module Psi
  class Deriv2
    # Mixin the InputGenerator
    include InputGenerator
    include Executor

    def initialize(task_obj)
      # Save which task we are created by
      set_task task_obj
      # Override the psi_module_name
      set_psi_module_name "Cints"
      # Set the generic command for cints
      set_binary_command Psi::Commands::DERIV2
    end
  end
  
  # Add cints ability to the Task class
  class Task
    def deriv2(*args)
      # convert to a hash
      args_hash = args[0]

      # Create a new input object
      input_obj = Psi::Deriv2.new self

      # Form the input hash and generate the input file
      input_hash = { "wfn" => wavefunction, "reference" => reference }
      input_hash = input_hash.merge(args_hash) unless args_hash == nil

      # Run the input module, sending the input file as keyboard input
      puts "2nd derivatives"
      input_obj.execute(input_hash)
    end
  end
end

# Create some global functions
# User can send additional input parameters to the function
def deriv2(*args)
  args_hash = args[0]
  Psi::global_task.deriv(args_hash)
end

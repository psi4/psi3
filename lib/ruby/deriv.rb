# Handle access to Psi cints module

module Psi
  class Deriv
    # Mixin the InputGenerator
    include InputGenerator
    include Executor

    def initialize(task_obj)
      # Save which task we are created by
      set_task task_obj
      # Override the psi_module_name
      set_psi_module_name "Cints"
      # Set the generic command for cints
      set_binary_command Psi::Commands::DERIV
    end
  end
  
  # Add cints ability to the Task class
  class Task
    def deriv(*args)
      # convert to a hash
      args_hash = args[0]

      # Create a new input object
      input_obj = Psi::Deriv.new self

      # Form the input hash and generate the input file
      input_hash = { "wfn" => wavefunction, "reference" => reference }
      input_hash = input_hash.merge(args_hash) unless args_hash == nil

      # Run the input module, sending the input file as keyboard input
      puts "derivatives"
      input_obj.execute(input_hash)
    end
  end
end

# Create some global functions
# User can send additional input parameters to the function
def deriv(*args)
  args_hash = args[0]
  Psi::global_task.deriv(args_hash)
end

# Handle access to Psi cints module

module Psi
  class Cints
    # Mixin the InputGenerator
    include InputGenerator
    include Executor

    def initialize(task_obj)
      # Save which task we are created by
      @task = task_obj
      # Set the generic command for cints
      set_binary_command Psi::Commands::INTS
    end
  end
  
  # Add cints ability to the Task class
  class Task
    def ints(*args)
      # convert to a hash
      args_hash = args[0]

      # Create a new input object
      input_obj = Psi::Cints.new self

      # Form the input hash and generate the input file
      input_hash = { }
      input_hash = input_hash.merge(args_hash) unless args_hash == nil

      # Run the input module, sending the input file as keyboard input
      puts "ints"
      input_obj.execute(input_hash)
    end

    def ints_derivative(*args)
      # convert to a hash
      args_hash = args[0]

      # Create a new input object
      input_obj = Psi::Cints.new self

      # Form the input hash and generate the input file
      input_hash = { }
      input_hash = input_hash.merge(args_hash) unless args_hash == nil

      # Run the input module, sending the input file as keyboard input
      puts "ints (1st derivative)"
      input_obj.execute(input_hash, Psi::Commands::DERIV)
    end
  end
end

# Create some global functions
# User can send additional input parameters to the function
def ints(*args)
  puts "in global ints"
  args_hash = args[0]
  Psi::global_task.ints(args_hash)
end

def ints_derivative(*args)
  args_hash = args[0]
  Psi::global_task.ints_derivative(args_hash)
end

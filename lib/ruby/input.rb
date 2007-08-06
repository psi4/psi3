# Handle access to Psi input module

module Psi
  class Input
    # Mixin the InputGenerator
    include InputGenerator
    include Executor
    
    # Static function used by all instances of Input
    def self.first_time_run
      if @input_first_time_run == nil
        @input_first_time_run = false
        Psi::Commands::INPUT
      else
        Psi::Commands::INPUTKEEP
      end
    end
    
    def initialize(task_obj)
      # Save which task we are created by
      @task = task_obj
      # Set the generic command
      set_binary_command Input.first_time_run
    end
  end
  
  # Add input ability to the Task class
  class Task
    # User can send additional input parameters to the function
    def input(*args)
      # convert to a hash
      args_hash = args[0]

      # Create a new input object
      input_obj = Psi::Input.new self

      # Make sure some label is set
      @label = "Default PSIRB label" if @label == nil

      # Form the input hash and generate the input file
      input_hash = { "label" => "\\\"#{@label.to_str}\\\"", "basis" => @basis }

      # Handle the geometry
      if @zmat != nil and @geometry == nil
        input_hash["zmat"] = @zmat
      elsif @zmat == nil and geometry != nil
        input_hash["geometry"] = @geometry
      elsif @zmat == nil and @geometry == nil
        puts "Error: Neither geometry nor zmat are set."
        exit 1
      else
        puts "Error: Both geometry and zmat are set. One must be nil."
        exit 1
      end
      input_hash = input_hash.merge(args_hash) unless args_hash == nil

      # Run the input module, sending the input file as keyboard input
      puts "input"
      input_obj.execute(input_hash)
    end
  end
end

# A generic global input function that call input on the global Task
def input(*args)
  args_hash = args[0]
  Psi::global_task.input(args_hash)
end
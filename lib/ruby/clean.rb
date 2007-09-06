# Handle access to Psi scf module

module Psi
  class Clean
    
    # Mixin the InputGenerator
    include InputGenerator
    include Executor
    
    def initialize(task_obj)
      @task = task_obj
      # Set the generic command for this class
      set_binary_command Psi::Commands::DONE
    end
  end
  
  # Add psiclean ability to Task
  class Task
    def clean(*args)
      # convert to a hash
      args_hash = args[0]

      # Create a new scf object
      clean_obj = Psi::Clean.new self

      # Form the input hash and generate the input file
      input_hash = { }

      # Merge what we've done with what the user wants
      input_hash = input_hash.merge(args_hash) unless args_hash == nil

      # Run the psiclean module, sending the input file as keyboard input
      puts "cleaning"
      clean_obj.execute(input_hash)
      
      # Return true
      true
    end
    
    def veryclean(*args)
      args_hash = args[0]
      
      # Run clean
      clean(args_hash)
      
      # Remove files in the local directory that match prefix.[0-9]*
      list = Dir[prefix + ".[0-9]*"]
      printf "Deleting: "
      list.each do |item|
        printf item + " "
        File.delete(item)
      end
      puts

      true
    end
  end
end

# Create some global functions
# User can send additional input parameters to the function
def clean(*args)
  # convert to a hash
  args_hash = args[0]
  Psi::global_task.clean(args_hash)
end

def veryclean(*args)
  # convert to a hash
  args_hash = args[0]
  Psi::global_task.veryclean(args_hash)
end

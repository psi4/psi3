# Handle access to Psi MP2 module

module Psi
  class MP2
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
      set_binary_command Psi::Commands::MP2
    end
  end
  
  # Add mp2 ability to Task
  class Task
    def mp2(*args)
      # convert to a hash
      args_hash = args[0]

      # Call transqt and ccsort unless told not to
      if args_hash != nil
        if args_hash.has_key?(:transqt) == false or 
           (args_hash.has_key?(:transqt) and args_hash[:transqt] == true)
          transqt(args_hash)
        end

        if args_hash.has_key?(:ccsort) == false or
           (args_hash.has_key?(:ccsort) and args_hash[:ccsort] == true)
          ccsort(args_hash)
        end
      else
        transqt(args_hash)
        ccsort(args_hash)
      end
      args_hash.delete(:transqt) unless args_hash == nil
      args_hash.delete(:ccsort)  unless args_hash == nil

      # Create a new mp2 object
      mp2_obj = Psi::MP2.new self

      # Form the input hash and generate the input file
      input_hash = { }

      # Check to see if the function arguments have the reference, if so use it, otherwise use
      # global setting
      if args_hash == nil or args_hash.has_key?("reference") == false
        input_hash["reference"] = reference
      end

      # Check the wavefunction
      if args_hash == nil or args_hash.has_key?("wfn") == false
        input_hash["wfn"] = wavefunction
      end

      # Merge what we've done with what the user wants
      input_hash = input_hash.merge(args_hash) unless args_hash == nil

      # Run the scf module, sending the input file as keyboard input
      puts "mp2"
      mp2_obj.execute(input_hash)
    end
  end
end

# Create some global functions
# User can send additional input parameters to the function
def mp2(*args)
  # convert to a hash
  args_hash = args[0]
  Psi::global_task.mp2(args_hash)
end

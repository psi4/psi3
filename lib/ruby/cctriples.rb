# Handle access to Psi scf module

module Psi
  class CCTriples
    
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
      set_binary_command Psi::Commands::CCTRIPLES
    end
  end
  
  # Add cctriples ability to Task
  class Task
    def ccsd_t(*args)
      # convert to a hash
      args_hash = args[0]

      # Ensure that ccsd is called unless told otherwise
      if args_hash != nil
        if args_hash.has_key?(:ccsd) == false or 
          (args_hash.has_key?(:ccsd) == true and args_hash[:ccsd] != false)
          ccsd(args_hash)
        end
      else
        ccsd(args_hash)
      end
      args_hash.delete(:ccsd)    unless args_hash == nil
      args_hash.delete(:transqt) unless args_hash == nil
      args_hash.delete(:ccsort)  unless args_hash == nil

      # Create a new scf object
      ccsd_t_obj = Psi::CCTriples.new self

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
      puts "ccsd(t)"
      ccsd_t_obj.execute(input_hash)
    end
  end
end

# Create some global functions
# User can send additional input parameters to the function
def ccsd_t(*args)
  # convert to a hash
  args_hash = args[0]
  Psi::global_task.ccsd_t(args_hash)
end

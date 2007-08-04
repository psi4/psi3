# Handle access to Psi scf module

module Psi
  class SCF
    
    # Mixin the InputGenerator
    include InputGenerator
    include Executor
    include SupportsAnalyticalGradients
    
    def initialize(task_obj)
      @task = task_obj
      # Set the generic command for this class
      set_binary_command Psi::Commands::SCF
    end
    
    def analytical_gradients(*args)
      # Need to call cints --deriv1
      ints_derivative(args[0])
    end
  end
  
  # Add scf ability to Task
  class Task
    def scf(*args)
      # convert to a hash
      args_hash = args[0]

      # Make sure cints is called
      if args_hash != nil 
        if args_hash.has_key?(:ints) == false or 
           (args_hash.has_key?(:ints) and args_hash[:ints] == true)
          ints("wfn" => "scf")
        end
      else
        ints("wfn" => "scf")
      end

      # Create a new scf object
      scf_obj = Psi::SCF.new self

      # Form the input hash and generate the input file
      input_hash = { }

      # Check to see if the user overrode the wavefunction
      if args_hash == nil or args_hash.has_key?("wfn") == false
        input_hash["wfn"] = wavefunction
      end
      args_hash.delete("wfn") unless args_hash == nil

      # Check to see if the user gave the reference wavefunction. if it is given, do not use the
      # global setting
      if args_hash == nil or args_hash.has_key?("reference") == false
        input_hash["reference"] = reference
      end
      args_hash.delete("reference") unless args_hash == nil

      # If we are doing analytical gradients need cscf to know
      if analytical_gradients == true and scf_obj.supports_analytical_gradients == true
        input_hash["dertype"] = 1
      else
        input_hash["dertype"] = "none"
      end

      # Merge what we've done with what the user wants
      input_hash = input_hash.merge(args_hash) unless args_hash == nil

      # Run the scf module, sending the input file as keyboard input
      puts "scf"
      scf_obj.execute(input_hash)

      # Check to see if we are supposed to compute analytical gradients
      if analytical_gradients == true and scf_obj.supports_analytical_gradients == true
        scf_obj.analytical_gradients(input_hash)
      end
    end

    def rhf(*args)
      # Convert to hash
      args_hash = args[0]

      # Set the Task reference to RHF
      # Must use @ sign here to access the instance variable, if we did reference then only
      # local variable is set
      @reference = "rhf"

      # Call scf
      scf(args_hash)
    end

    def rohf(*args)
      args_hash = args[0]

      # Set the Task reference to ROHF
      @reference = "rohf"

      # Call scf
      scf(args_hash)
    end

    def uhf(*args)
      args_hash = args[0]

      # Set the global reference to UHF
      @reference = "uhf"

      # Call scf
      scf(args_hash)
    end    
  end
end

# Create some global functions
# User can send additional input parameters to the function
def scf(*args)
  # convert to a hash
  args_hash = args[0]
  Psi::global_task.scf(args_hash)
end

def rhf(*args)
  # Convert to hash
  args_hash = args[0]
  Psi::global_task.rhf(args_hash)
end

def rohf(*args)
  args_hash = args[0]
  Psi::global_task.rohf(args_hash)
end

def uhf(*args)
  args_hash = args[0]
  Psi::global_tast.uhf(args_hash)
end

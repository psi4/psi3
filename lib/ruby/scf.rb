# Handle access to Psi scf module
require 'pp'

module Psi
  class SCF
    include InputGenerator
    include Executor
    # Which references does this method support
    def self.supports_analytic_gradients
      { "rhf" => true, "rohf" => true, "uhf" => true, "twocon" => true }
    end
    def self.supports_analytic_second_derivatives
      { "rhf" => true, "rohf" => false, "uhf" => false, "twocon" => false }
    end
    
    def initialize(task_obj)
      @task = task_obj
      # Set the generic command for this class
      set_binary_command Psi::Commands::SCF
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
      else
        input_hash["wfn"] = args_hash["wfn"]
      end
      args_hash.delete("wfn") unless args_hash == nil

      # Check to see if the user gave the reference wavefunction. if it is given, do not use the
      # global setting
      if args_hash == nil or args_hash.has_key?("reference") == false
        input_hash["reference"] = reference
      else
        input_hash["reference"] = args_hash["reference"]
      end
      args_hash.delete("reference") unless args_hash == nil

      # If we are doing analytical gradients then cscf needs to know
      if get_gradients == true and SCF.supports_analytic_gradients[reference] == true
        input_hash["dertype"] = "first"
      elsif get_second_deriv == true and SCF.supports_analytic_second_derivatives[reference] == true
        input_hash["dertype"] = "second"
      else
        input_hash["dertype"] = "none"
      end

      # Merge what we've done with what the user wants
      input_hash = input_hash.merge(args_hash) unless args_hash == nil

      # Run the scf module, sending the input file as keyboard input
      puts "scf"
      scf_obj.execute(input_hash)

      # Check to see if we are supposed to compute analytical gradients, only do if
      #  1. It is supported
      #  2. The user requested it
      #  3. The wavefunction is SCF
      if get_gradients == true and SCF.supports_analytic_gradients[reference] == true and wavefunction.casecmp("scf") == 0
        # Make sure the wavefunction is set
        input_hash["wfn"] = "scf"
        deriv(input_hash)
      elsif get_second_deriv == true and SCF.supports_analytic_second_derivatives[reference] == true and wavefunction.casecmp("scf") == 0
        # Make sure the wavefunction is set
        input_hash["wfn"] = "scf"
        transqt(input_hash)
        deriv2(input_hash)
        propint(input_hash)
        cphf(input_hash)
      end
      
      # Return the total energy to the user
      Psi::Chkpt::etot
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
  Psi::global_task.uhf(args_hash)
end

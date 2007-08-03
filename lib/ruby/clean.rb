# Handle access to Psi scf module

module Psi
  class Clean
    
    # Mixin the InputGenerator
    include InputGenerator
    include Executor
    include NoAnalyticalGradients
    
    def initialize
      # Set the generic command for this class
      set_binary_command Psi::Commands::DONE
    end
    
    def analytical_gradients(*args)
      # Need to call cints --deriv1
      ints_derivative(args[0])
    end
  end
end

# Create some global functions
# User can send additional input parameters to the function
def clean(*args)
  # convert to a hash
  args_hash = args[0]
    
  # Create a new scf object
  clean_obj = Psi::Clean.new
  
  # Form the input hash and generate the input file
  input_hash = { }
    
  # Merge what we've done with what the user wants
  input_hash = input_hash.merge(args_hash) unless args_hash == nil
  
  # Run the scf module, sending the input file as keyboard input
  puts "cleaning"
  clean_obj.execute(input_hash)
end

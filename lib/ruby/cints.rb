# Handle access to Psi cints module

module Psi
  class Cints
    # Mixin the InputGenerator
    include InputGenerator
    include Executor

    def initialize
      # Set the generic command for cints
      set_binary_command Psi::Commands::INTS
    end
  end
end

# Create some global functions
# User can send additional input parameters to the function
def ints(*args)
  # convert to a hash
  args_hash = args[0]
  
  # Create a new input object
  input_obj = Psi::Cints.new
  
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
  input_obj = Psi::Cints.new
  
  # Form the input hash and generate the input file
  input_hash = { }
  input_hash = input_hash.merge(args_hash) unless args_hash == nil
    
  # Run the input module, sending the input file as keyboard input
  puts "ints (1st derivative)"
  input_obj.execute(input_hash, Psi::Commands::DERIV)
end

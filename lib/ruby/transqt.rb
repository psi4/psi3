# Handle access to Psi scf module

module Psi
  class TransQT
    
    # Mixin the InputGenerator
    include InputGenerator
    include Executor
    
    def initialize
      # Set the generic command for this class
      set_binary_command Psi::Commands::TRANSQT
    end
  end
end

# Create some global functions
# User can send additional input parameters to the function
def transqt(*args)
  # convert to a hash
  args_hash = args[0]
    
  # Create a new scf object
  transqt_obj = Psi::TransQT.new
  
  # Form the input hash and generate the input file
  input_hash = { }
  
  # Check to see if the function arguments have the reference, if so use it, otherwise use
  # global setting
  if args_hash == nil or args_hash.has_key?("reference") == false
    input_hash["reference"] = Psi::reference
  end
  
  # Check the wavefunction
  if args_hash == nil or args_hash.has_key?("wfn") == false
    input_hash["wfn"] = Psi::wavefunction
  end
  
  # Are we doing a back transformation?
  if args_hash != nil and args_hash.has_key?(:backtr) and args_hash[:backtr] == true
    transqt_obj.set_binary_command Psi::Commands::BACKTRANSQT
  end
  args_hash.delete(:backtr) unless args_hash == nil
  
  # Merge what we've done with what the user wants
  input_hash = input_hash.merge(args_hash) unless args_hash == nil
  
  # Run the scf module, sending the input file as keyboard input
  puts "transqt"
  transqt_obj.execute(input_hash)
end

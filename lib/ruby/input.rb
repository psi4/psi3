# Handle access to Psi input module

module Psi
  class Input
    # Mixin the InputGenerator
    include InputGenerator
    include Executor
    
    def self.first_time_run
      if @input_first_time_run == nil
        @input_first_time_run = false
        Psi::Commands::INPUT
      else
        Psi::Commands::INPUTKEEP
      end
    end
    
    def initialize
      # Set the generic command
      set_binary_command Input.first_time_run
    end
    
  end
end

# Create some global functions
# User can send additional input parameters to the function
def input(*args)
  # convert to a hash
  args_hash = args[0]
  
  # Create a new input object
  input_obj = Psi::Input.new
  
  # Make sure some label is set
  Psi::label = "Default PSIRB label" if Psi::label == nil
  
  # Form the input hash and generate the input file
  input_hash = { "label" => "\\\"#{Psi::label.to_str}\\\"", "basis" => Psi::basis }
  
  # Handle the geometry
  geometry = Psi::geometry
  zmat = Psi::zmat
  if zmat != nil and geometry == nil
    input_hash["zmat"] = zmat
  elsif zmat == nil and geometry != nil
    input_hash["geometry"] = geometry
  elsif zmat == nil and geometry == nil
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

# Handle access to optking module

module Psi
  class OptKing
    
    PRE  = 1
    STEP = 2
    POST = 3
    
    # Mixin the InputGenerator
    include InputGenerator
    # We don't mix in the Executor because optking returns error codes corresponding to certain actions
    
    def execute(input_file, program_name)
      puts "DEBUG: For #{binary} would have used:\n#{input_file.string}" if Psi::check_commands == true
      
      if Psi::check_commands != true
        # Run optking with our input file
        `echo "#{input_file.string}" | #{binary} -f -`
        status = $?
        
        # Return the status code to the caller so that it can handle what to do next
        return status
      end
    end
  end
end

# User is expected to provide a code block that tells optking what to do
# for a single step.
def optking(*args)
  
  # This code is not complete. Error out
  puts "Error: The optking interface does not work."
  exit 1
  
  # Convert to a hash
  args_hash = args[0]
  
  # Create a new optking object
  optking_obj = Psi::OptKing.new
  
  # Form the input hash and generate input file
  input_hash = {}
  input_hash = input_hash.merge(args_hash) unless args_hash == nil
  
  # Generate input file
  input_file = optking_obj.generate_input input_hash

  puts "Optimizing geometry..."
  Psi::indent_puts

  # Step 0, run input to generate initial checkpoint file contents
  input(args_hash)
  
  # Step 1, run a single point so that optking has something to work with
  result = yield Psi::OptKing::PRE
  Psi::unindent_puts
  
  if result != nil and result != 0
    puts "Problem occurred during initial single point."
    exit 1
  end
end
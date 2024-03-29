# Code that simplifies computing a potential energy surface

# Add pes ability to Task
module Psi
  class Task
    # Function arguments:
    #  :num_disps => #          (required)
    #  :report => true or false (default => true)
    #  :col_labels => array     (required, column headings)
    #  :file => string          (optional)
    #  :format => formatstring  (optional) if "csv" saves the file in csv format.
    #
    # Expects a yield block that runs the single point and returns the information to collect.
    #
    # Returns the information that was collected
    def pes(*args)
      # Convert the array to a hash
      args_hash = args[0]

      format = "normal"
      format = args_hash[:format] if args_hash.has_key?(:format)
      
      # Make sure num_disps is provided
      if args_hash.has_key?(:num_disps) == false
        puts "Error: In function \"pes\", :num_disps required."
        exit 1
      end

      # Check to see if :report is given and if it is true :col_label is required
      if args_hash.has_key?(:report) and args_hash[:report] == true
        if args_hash.has_key?(:col_labels) == false
          puts "Error: In function \"pes\", :col_labels is required if :report => true."
          exit 1
        end
      end

      # Create the result array
      result = []

      # Loop over the displacements
      puts "Beginning PES scan..."
      Psi::indent_puts
      (1..args_hash[:num_disps]).each do |ndisp|
        puts "Displacement #{ndisp}"
        Psi::indent_puts
        result += yield(ndisp)
        Psi::unindent_puts
      end
      Psi::unindent_puts

      # Report the results
      if args_hash.has_key?(:report) and args_hash[:report] == true
        puts "Results from PES scan:"
        args_hash[:col_labels].each do |label|
          printf(sprintf(" %-14s ", label.to_s))
        end
        printf "\n"

        col_count = (args_hash[:col_labels].length) - 1
        row_count = (result.length / args_hash[:col_labels].length) - 1

        (0..row_count).each do |row|
          (0..col_count).each do |col|
            indx = row * (col_count+1) + col
            if result[indx].instance_of? Fixnum or result[indx].instance_of? Bignum
              printf(sprintf(" %-14d ", result[indx]))
            elsif result[indx].instance_of? Float
              printf(sprintf(" %-14.9f ", result[indx]))
            else
              printf(sprintf(" %-14s ", result[indx].to_s))
            end
          end
          printf "\n"
        end
      end

      # Does the user want the data saved to a file?
      if args_hash.has_key?(:file)
        File.open(args_hash[:file], "w") do |file|
          args_hash[:col_labels].each do |label|
            if format == "normal"
              file.printf(sprintf(" %-14s ", label.to_s))
            elsif format == "csv"
              file.printf(label.to_s + ",")
            end
          end
          file.printf "\n"

          col_count = (args_hash[:col_labels].length) - 1
          row_count = (result.length / args_hash[:col_labels].length) - 1

          (0..row_count).each do |row|
            (0..col_count).each do |col|
              indx = row * (col_count+1) + col
              if format == "normal"
                if result[indx].instance_of? Fixnum or result[indx].instance_of? Bignum
                  file.printf(sprintf(" %-14d ", result[indx]))
                elsif result[indx].instance_of? Float
                  file.printf(sprintf(" %-14.9f ", result[indx]))
                else
                  file.printf(sprintf(" %-14s ", result[indx].to_s))
                end
              elsif format == "csv"
                file.printf(result[indx].to_s + ",")
              end
            end
            file.printf "\n"
          end
        end
      end

      # If the user provided :col_labels, convert result into a more 
      # manageable form
      if args_hash.has_key?(:col_labels)
        new_result = {}
        # Add empty arrays to the new hash
        args_hash[:col_labels].each do |label|
          new_result[label] = []
        end

        # Go through added the data to the appropriate array
        col_count = (args_hash[:col_labels].length) - 1
        row_count = (result.length / args_hash[:col_labels].length) - 1

        (0..row_count).each do |row|
          (0..col_count).each do |col|
            new_result[args_hash[:col_labels][col]] << result[row * (col_count+1) + col]
          end
        end

        result = new_result
      end

      return result
    end
  end
end

# Provide global ability to pes
def pes(*args)
  args_hash = args[0]
  Psi::global_task.pes(args_hash) do |disp|
    yield disp
  end
end

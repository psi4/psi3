# Global include file for Psi3
#  Includes the basics to get you started

# Other Ruby modules that are required
require 'stringio'

# Please note that not all the objects are given here. Any commands that communicate directly to
# any Psi3 library function is implemented in the source for psirb.
module Psi
  
  module Commands
    INPUT       = "input"
    INPUTKEEP   = "input --keepoutput"
    UINPUT      = "input --chkptgeom"
    INIT        = INPUT
    DONE        = "psiclean"
    INTS        = "cints"
    SCF         = "cscf"
    LOCALIZE    = "localize"
    DERIV       = "cints --deriv1"
    DERIV2      = "cints --oeprop"
    PROPINT     = "cins --oeprop"
    TRANSQT     = "transqt"
    CCTRANS     = "transqt2"
    BACKTRANSQT = "transqt --backtr"
    CPHF        = "cphf"
    RESPONSE    = "response"
    NORMCO      = "normco"
    OEPROP      = "oeprop"
    OPTKING     = "optking"
    STABLE      = "stable"
    CIS         = "cis"
    UGEOM       = "ugeom"
    CCSORT      = "ccsort"
    CCRESET     = "ccsort --reset"
    CCENERGY    = "ccenergy"
    CCTRIPLES   = "cctriples"
    CCHBAR      = "cchbar"
    CCLAMBDA    = "cclambda"
    CCEOM       = "cceom"
    CCRESPONSE  = "ccresponse"
    CCDENSITY   = "ccdensity"
  end
  
  def self.check_commands=(val)
    @check_commands = val
  end
  
  def self.check_commands
    @check_commands
  end
  
  # Routines to set active geometry, referred to only by Input
  def self.geometry=(val)
    @geometry = val
  end
  def self.geometry
    @geometry
  end
  
  # Routines to set active zmat, referred to only by Input
  def self.zmat=(val)
    $zmat = val
  end
  def self.zmat
    $zmat
  end
  
  # Routines to remember basis
  def self.basis=(val)
    @basis = val
  end
  def self.basis
    @basis
  end
  
  # Set the label
  def self.label=(val)
    @label = val
  end
  def self.label
    @label
  end
  
  # Set the label
  def self.wavefunction=(val)
    @wavefunction = val
  end
  def self.wavefunction
    if @wavefunction == nil
      # Was never set, assume SCF
      puts "Warning: Wavefunction not set. Assuming SCF."
      return "SCF"
    end
    @wavefunction
  end
  
  # Set the flag for analytical gradients
  def self.analytical_gradients=(val)
    @analytical_gradients = val
  end
  def self.analytical_gradients
    if @analytical_gradients == nil
      return false
    else
      @analytical_gradients
    end
  end
  
  # Set the global reference
  def self.reference=(val)
    @reference = val
  end
  def self.reference
    if @reference == nil
      # We don't know what the user wanted to do.
      # Tell the user about it:
      puts "Warning: Reference not set. Assuming RHF."
      return "RHF"
    end
    @reference
  end
  
  module InputGenerator
    def generate_value(item)
      result = StringIO.new("", "w")
      if item.kind_of?(Array)
        # Encompass an array with ()
        result.printf "("
        item.each do |x|
          sub_result = generate_value(x)
          result.printf(sub_result.string)
        end
        # close off the array
        result.puts ")"
      elsif item.kind_of?(Hash)
        item.each do |key,value|
          result.puts "#{key} = #{generate_value(value).string}"
        end
      else
        # Tell everything else to convert to a string
        result.printf "#{item.to_s} "
      end
      
      return result
    end
    
    #
    # Set the input parameters that are needed to run the given psi module
    # Expects a hash containing properties to be set
    #  { "jobtype" => "sp", "dertype" => "none" }
    # Uses the name of the encompassing class to name the module for psi3 input file
    def generate_input(input_hash)
      input = StringIO.new("", "w")
      psi_module_name = self.class.name.slice!(5..-1)
      
      # Make sure the scratch location stuff is included
      input.puts "psi:("
      input.puts "  files:("
      input.puts "    default:("
      input.puts "      nvolume=1"
      input.puts "      volume1=\\\"#{Psi::scratch}\\\""
      input.puts "    )"
      input.puts "    file32: (nvolume=1 volume1=\\\"./\\\")"
      input.puts "  )"
      input.puts ")"
      
      # Begin the input
      input.puts "#{psi_module_name}:("
      
      # Go through the hash and add entries into the input
      sub_result = generate_value(input_hash)
      input.puts(sub_result.string)
      # Close off the input
      input.puts ")"
      
      return input
    end
  end
  
  module Executor
    def set_binary_command(binary)
      @program = binary
    end
    
    def get_binary_command
      @program
    end
    
    def execute_internal(input_file, binary=nil)
      puts "DEBUG: For #{binary} would have used:\n#{input_file.string}" if Psi::check_commands == true and binary != nil
      if Psi::check_commands == true and get_binary_command != nil
        puts "DEBUG: For #{get_binary_command} would have used:\n#{input_file.string}"
      end
      
      if Psi::check_commands != true
        prefix = ""
        prefix = "-p #{Psi::prefix}" if Psi::prefix != nil
        
        if binary != nil
#         `echo "#{input_file.string}" | #{binary} -f - >& /dev/null`
          `echo "#{input_file.string}" | #{binary} #{prefix} -f -`
        elsif get_binary_command != nil
#         `echo "#{input_file.string}" | #{get_binary_command} -f - >& /dev/null`
          `echo "#{input_file.string}" | #{get_binary_command} #{prefix} -f -`
        else
          # TODO: Add coloring commands to text
          puts "Error: Executor.execute: Unable to determine which module to run."
          exit 1
        end
        # Check to see if the above line worked
        status = $?
#       puts "Program executed successfully." if status == 0
#       puts "Program executed un-successfully." if status != 0
      
        # Note for optking a specialized execute is needed
        if status != 0
          puts "Module exited with an error."
          puts "Used:\n#{input_file.string}"
          exit 1
        end
      end
    end
    
    # This function requires that InputGenerator be used
    def execute(input_hash, binary=nil)
      # Form the input file
      input_file = generate_input input_hash

      # execute the module
      execute_internal(input_file, binary)
    end
  end
  
  module SupportsAnalyticalGradients
    def supports_analytical_gradients
      true
    end
  end
  
  module NoAnalyticalGradients
    def supports_analytical_gradients
      false
    end
  end
end

require 'input'
require 'cints'
require 'scf'
require 'optking'
require 'clean'
require 'pes'
require 'transqt'
require 'ccsort'
require 'ccenergy'
require 'cctriples'
require 'symbols'
  
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

  # Access to the global task object. The first time it is retrieved the object is created.
  def self.global_task=(val)
    @global_task = val
  end
  def self.global_task
    if @global_task == nil
      # Creat the global task
      @global_task = Psi::Task.new
    end
    @global_task
  end

  def self.prefix=(val)
    Psi::global_task.prefix=val
  end
  def self.prefix
    Psi::global_task.prefix
  end
  
  def self.scratch=(val)
    Psi::global_task.scratch=val
  end
  def self.scratch
    Psi::global_task.scratch
  end
  
  def self.check_commands=(val)
    Psi::global_task.check_commands = val
  end
  
  def self.check_commands
    Psi::global_task.check_commands
  end
  
  # Routines to set active geometry, referred to only by Input
  def self.geometry=(val)
    Psi::global_task.geometry = val
  end
  def self.geometry
    Psi::global_task.geometry
  end
  
  # Routines to set active zmat, referred to only by Input
  def self.zmat=(val)
    Psi::global_task.zmat = val
  end
  def self.zmat
    Psi::global_task.zmat
  end
  
  # Routines to remember basis
  def self.basis=(val)
    Psi::global_task.basis = val
  end
  def self.basis
    Psi::global_task.basis
  end
  
  # Set the label
  def self.label=(val)
    Psi::global_task.label = val
  end
  def self.label
    Psi::global_task.label
  end
  
  # Set the label
  def self.wavefunction=(val)
    Psi::global_task.wavefunction = val
  end
  def self.wavefunction
    Psi::global_task.wavefunction
  end
  
  # Set the flag for analytical gradients
  def self.analytical_gradients=(val)
    Psi::global_task.analytical_gradients=val
  end
  def self.analytical_gradients
    Psi::global_task.analytical_gradients
  end
  
  # Set the global reference
  def self.reference=(val)
    Psi::global_task.reference = val
  end
  def self.reference
    Psi::global_task.reference
  end
  
  # Additions to the Task class
  class Task
    def check_commands=(val)
      @check_commands = val
    end
    def check_commands
      if @check_commands == nil
        @check_commands = false
      end
      @check_commands
    end
    
    def reference=(val)
      @reference = val
    end
    def reference
      if @reference == nil
        puts "Warning: Reference not set. Assuming RHF."
        return "RHF"
      end
      @reference
    end

    def wavefunction=(val)
      @wavefunction = val
    end
    def wavefunction
      if @wavefunction == nil
        # Was never set, assume SCF
        puts "Warning: Wavefunction not set. Assuming SCF."
        return "SCF"
      end
      @wavefunction
    end
    
    # Routines to set active geometry, referred to only by Input
    def geometry=(val)
      @geometry = val
    end
    def geometry
      @geometry
    end

    # Routines to set active zmat, referred to only by Input
    def zmat=(val)
      @zmat = val
    end
    def zmat
      @zmat
    end

    # Routines to remember basis
    def basis=(val)
      @basis = val
    end
    def basis
      @basis
    end

    # Set the label
    def label=(val)
      @label = val
    end
    def label
      @label
    end
    
    def analytical_gradients=(val)
      @analytical_gradients = val
    end
    def analytical_gradients
      if @analytical_gradients == nil
        return false
      end
      @analytical_gradients
    end
    
    def create(*args)
      args_hash = args[0]
      t = Psi::Task.new args_hash
      if t == nil
        puts "Error: Unable to create Task object."
        exit 1
      end
      if block_given?
        yield t
      else
        t
      end
    end
  end # of Task
  
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
      input.puts "      volume1=\\\"#{@task.scratch}\\\""
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
        prefix = "-p #{@task.prefix}" if @task.prefix != nil
        
        if binary != nil
#         `echo "#{input_file.string}" | #{binary} -f - >& /dev/null`
          `echo "#{input_file.string}" | #{binary} #{prefix} -f -`
        elsif get_binary_command != nil
#         `echo "#{input_file.string}" | #{get_binary_command} -f - >& /dev/null`
          `echo "#{input_file.string}" | #{get_binary_command} #{prefix} -f -`
        else
          puts "Error: Executor.execute: Unable to determine which module to run."
          exit 1
        end
        # Check to see if the above line worked
        status = $?
      
        # Note for optking a specialized execute is needed
        if status != 0
          puts "Module exited with an error."
          puts "Error code: #{status}"
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
require 'chkpt'

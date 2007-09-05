module Psi
  module Chkpt

  	#rb_define_module_function(rubyChkpt, "exist?",  RUBYCAST(ruby_psi_chkpt_exist),     1);
  	def self.exist?(key)
  	  Psi::global_task.exist?(key)
	  end
	  
  	#rb_define_module_function(rubyChkpt, "exists?", RUBYCAST(ruby_psi_chkpt_exist),     1);
  	def self.exists?(key)
  	  Psi::global_task.exists?(key)
	  end
	  
  	#rb_define_module_function(rubyChkpt, "label", 	RUBYCAST(ruby_psi_chkpt_label_get), 0);
  	def self.label
  	  Psi::global_task.label
	  end
	  
  	#rb_define_module_function(rubyChkpt, "escf",	RUBYCAST(ruby_psi_chkpt_escf_get),  0);
  	def self.escf
  	  Psi::global_task.escf
	  end
	  
  	#rb_define_module_function(rubyChkpt, "escf=",	RUBYCAST(ruby_psi_chkpt_escf_set),  1);
  	def self.escf=(val)
  	  Psi::global_task.escf=val
	  end
	  
  	#rb_define_module_function(rubyChkpt, "eref",	RUBYCAST(ruby_psi_chkpt_eref_get),  0);
  	def self.eref
  	  Psi::global_task.eref
	  end
	  
  	#rb_define_module_function(rubyChkpt, "eref=", 	RUBYCAST(ruby_psi_chkpt_eref_set),  1);
  	def self.eref=(val)
  	  Psi::global_task.eref=val
	  end
	  
  	#rb_define_module_function(rubyChkpt, "ecorr", 	RUBYCAST(ruby_psi_chkpt_ecorr_get), 0);
  	def self.ecorr
  	  Psi::global_task.ecorr
	  end
	  
  	#rb_define_module_function(rubyChkpt, "ecorr=", 	RUBYCAST(ruby_psi_chkpt_ecorr_set), 1);
  	def self.ecorr=(val)
  	  Psi::global_task.ecorr=val
	  end
  	
  	#rb_define_module_function(rubyChkpt, "enuc", 	RUBYCAST(ruby_psi_chkpt_enuc_get),  0);
  	def self.enuc
  	  Psi::global_task.enuc
	  end
	  
  	#rb_define_module_function(rubyChkpt, "enuc=", 	RUBYCAST(ruby_psi_chkpt_enuc_set),  1);
  	def self.enuc=(val)
  	  Psi::global_task.enuc=val
	  end
	  
  	#rb_define_module_function(rubyChkpt, "efzc", 	RUBYCAST(ruby_psi_chkpt_efzc_get),  0);
  	def self.efzc
  	  Psi::global_task.efzc
	  end
	  
  	#rb_define_module_function(rubyChkpt, "efzc=", 	RUBYCAST(ruby_psi_chkpt_efzc_set),  1);
  	def self.efzc=(val)
  	  Psi::global_task.efzc=val
	  end
	  
  	#rb_define_module_function(rubyChkpt, "etot", 	RUBYCAST(ruby_psi_chkpt_etot_get),  0);
  	def self.etot
  	  Psi::global_task.etot
	  end
	  
  	#rb_define_module_function(rubyChkpt, "etot=", 	RUBYCAST(ruby_psi_chkpt_etot_set),  1);
  	def self.etot=(val)
  	  Psi::global_task.etot=val
	  end
	  
  	#rb_define_module_function(rubyChkpt, "disp", 	RUBYCAST(ruby_psi_chkpt_etot_get),  0);
  	def self.disp
  	  Psi::global_task.disp
	  end
	  
  	#rb_define_module_function(rubyChkpt, "disp=", 	RUBYCAST(ruby_psi_chkpt_etot_set),  1);
  	def self.disp=(val)
  	  Psi::global_task.disp=val
	  end
	  
  	#rb_define_module_function(rubyChkpt, "eccsd",   RUBYCAST(ruby_psi_chkpt_eccsd_get), 0);
  	def self.eccsd
  	  Psi::global_task.eccsd
	  end
	  
  	#rb_define_module_function(rubyChkpt, "e_t",     RUBYCAST(ruby_psi_chkpt_e_t_get),   0);
  	def self.e_t
  	  Psi::global_task.e_t
	  end
	  
	  def self.emp2
	    Psi::global_task.emp2
    end
    
    def self.num_irreps
      Psi::global_task.num_irreps
    end
  end
end

CC = gcc -g  -m64 -pedantic
OPT = -O2 
ACNUC_CSRC = /panhome/banques/csrc/trunk
OBJS_DIR = /panhome/penel/penel/svn/csrc/pipelinemaj
headers = $(ACNUC_CSRC)/dir_acnuc.h
BIN_DIR = $(shell if [ `uname` = SunOS ]; then echo /panhome/penel/penel/svn/bin/SunOS; else echo /panhome/penel/penel/svn/bin/Debian ; fi)
CFLAGS  = $(OPT) -I$(ACNUC_CSRC)
LIBACNUC = $(shell if [ `uname` = SunOS ]; then echo cacnuc64; else echo cacnucdeb; fi)
LIBNAME = maj
LIBFULLNAME = lib$(LIBNAME).a


# --------------------------
# Definition des executables
# --------------------------
CPROGS = add_fam_annot_genomes
# dupl_entry_no_para  test_sup_new_mod_hov genere_short_cds_GR_ENS \
#	    count_masked_residues add_fam_annot_genomes add_fam_annot_sp\
#	    get_def_hov  convert_ali convert_phylip listmn_2_fasta\
#	    link_sp2embl reformat_sp2embl\
#	    add_fam_annot_embl generate_CLUSTAL_data get_options_gb\
#	    phyml2nhx run_addroot annot_tree

# ---------------------
# Definition des objets
# ---------------------
LIB_OBJS =  $(OBJS_DIR)/misc.o $(OBJS_DIR)/convert_nexus.o $(OBJS_DIR)/getnseq.o

# --------------------------
# --------------------------
# Rules
# --------------------------
# --------------------------

# Dans le all: mettre la liste des programes en a compiler en 32 et $(CALLPROGS64) 
# =============================================================================

all :   $(CPROGS)
	@echo " "
	@echo "******************************************"
	@echo "Compilation de toutes les sources terminée"
	@echo "******************************************"

	@echo " "



$(LIB_OBJS) : $(headers) 	


$(LIBFULLNAME) : $(LIB_OBJS)
	ar rv $(LIBFULLNAME) $(LIB_OBJS)
	
	
# Compilation des executables
# ----------------------------

#$(CPROGS) : $(LIBFULLNAME) ; @echo $@.c
$(CPROGS) : $(LIBFULLNAME) barr
	@echo " "
	@echo "Compilation de $@  avec $(LIBACNUC)"
	$(CC)   -o $@ $@.c -L. -l$(LIBNAME)  -L$(ACNUC_CSRC) -l$(LIBACNUC) -lm
#	@echo "Copie $@ dans $(BIN_DIR)"	
#	cp $@ $(BIN_DIR)/$@
barr : ; @echo $@.c


## Implicit actions
.c.o :	
	@echo " "
	@echo "Compilation implicite"
	$(CC)  -c $(CFLAGS) $<
	

clean:
	@echo "Nettoyage des fichiers intermediaires"
	rm -rf *.o
	

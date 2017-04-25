/* 

###

Copyright CNRS 2008. Contributor: Simon Penel simon.penel@univ-lyon1.fr

###

This software is a computer program whose purpose is to annotate protein/protein-coding gene sequences
within an ACNUC sequence database with homologous gene family identifiers, as would have been generated 
by a sequence homology clustering algortihm. This is a key part of production of HOGENOM databases 
(Penel et al. 2009 BMC Bioinformatics, 10(S6):S3).

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.

####

Ce logiciel est un programme informatique servant à [rappeler les
caractéristiques techniques de votre logiciel]. 

Ce logiciel est régi par la licence CeCILL soumise au droit français et
respectant les principes de diffusion des logiciels libres. Vous pouvez
utiliser, modifier et/ou redistribuer ce programme sous les conditions
de la licence CeCILL telle que diffusée par le CEA, le CNRS et l'INRIA 
sur le site "http://www.cecill.info".

En contrepartie de l'accessibilité au code source et des droits de copie,
de modification et de redistribution accordés par cette licence, il n'est
offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,
seule une responsabilité restreinte pèse sur l'auteur du programme,  le
titulaire des droits patrimoniaux et les concédants successifs.

A cet égard  l'attention de l'utilisateur est attirée sur les risques
associés au chargement,  à l'utilisation,  à la modification et/ou au
développement et à la reproduction du logiciel par l'utilisateur étant 
donné sa spécificité de logiciel libre, qui peut le rendre complexe à 
manipuler et qui le réserve donc à des développeurs et des professionnels
avertis possédant  des  connaissances  informatiques approfondies.  Les
utilisateurs sont donc invités à charger  et  tester  l'adéquation  du
logiciel à leurs besoins dans des conditions permettant d'assurer la
sécurité de leurs systèmes et ou de leurs données et, plus généralement, 
à l'utiliser et l'exploiter dans les mêmes conditions de sécurité. 

Le fait que vous puissiez accéder à cet en-tête signifie que vous avez 
pris connaissance de la licence CeCILL, et que vous en avez accepté les
termes.

####

Requires dir_acnuc header, to be found at: http://doua.prabi.fr/databases/acnuc_data/dir_acnuc
Interoperates with ACNUC database, see http://doua.prabi.fr/databases/acnuc

Annote les cds dans les génomes.
--------------------------------

a faire ; nb de familles =106 au liu de 105...
cc -o add_fam_annot_ensembl add_fam_annot_ensembl.c -I/pbil2/banques/csrc/ -L/pbil2/banques/csrc/ -lcacnucsol
ou utiliser le Makefile


Septembre 2007 :

Je ne compte plus les bits, mais j'utilise simplement la valeur de faddr  garace  next_annots64(&faddr)

Bug decouvert:
quand un CDS est decrit sur plusieurs lignes ca plante:

exeemple DANRE21_54
on obtient:

FT   CDS             join(complement(28005..28051),complement(27411..27787),
FT                   /gene_family="HBG071438"
FT                   complement(22985..23107),complement(22731..22807),
FT                   complement(20873..20973),complement(18698..18851),
FT                   complement(7110..7259),complement(6325..6390),
FT                   complement(4737..4838),complement(3144..3278))
ce qui fait que le cds n,'est pas construit dans acnuc.
Il faut ecrire gene family appres la dernier ligne complement

tester :
FT                   complement(2495..3278))
et
FT                   (



Autre probleme:
FT                   /temporary_systematic_id="Tb11.18.0008"
FT   misc_feature    join(602757..602825,602853..602921,602976..603044,603072..603140,603255..603323,603336..603404,603465..603524,603552..603611,603669..603737,603780..603848,603909..603968,603996..604064,604212..604280,604584..604652,604878..604937,6049
65..605033,605094..605162,605175..605243,605304..605363,605391..605459,605760..605864,605907..605975,606012..606080,606108..606176,606195..606263)
FT                   /algorithm="TMHMM 2.0"
FT                   /colour=0

est traduit par 

FT                   /temporary_systematic_id="Tb11.18.0008"
FT   misc_feature    join(602757..602825,602853..602921,602976..603044,603072..603140,603255..603323,603336..603404,603465..603524,603552..603611,603669..603737,603780..603848,603909..603968,603996..604064,604212..604280,604584..604652,604878..604937,6049

65..605033,605094..605162,605175..605243,605304..605363,605391..605459,605760..605864,605907..605975,606012..606080,606108..606176,606195..606263)
FT                   /algorithm="TMHMM 2.0"
FT                   /colour=0



*/

#include "/panhome/banques/csrc/trunk/dir_acnuc.h"




struct SEQ {
	char cds[WIDTH_MAX];	/*Nom du CDS*/
	char mere[WIDTH_MAX];	/*Nom de la mere*/
	char fam[WIDTH_MAX];	/*Nom de la famille*/	
	off_t adresse;		/*adresse dans le fichier plat*/
	off_t madresse;		/*adresse dans le fichier plat par rapport a la mere*/
	} *cdsfam;
int nb_seq_tot;

int nb_warning;

void load_genefam(FILE *fic);
char *check_alloc(int nbrelt, int sizelt);
int sort_cdsfam(const void *a ,const void *b);
int add_annot(char *gene, FILE *outfile, int *nbfam);
int print_fam(int cur_cds, FILE *outfile);
void do_traitement(FILE *mn_file, FILE *outfile);



/************************* print_mess *************************************/
void print_mess(void)
{
fprintf(stderr, 
"                                                                        \n\
add_fam_annot_genomes: reads the  EMBL ACNUC database, and create           \n\
a new flat file with the gene_family identifiers in the CDS              \n\
annotations.                                                             \n\
                                                                         \n\
Usage: add_fam_annot_genomes list_mn GENE_FAM_file outfile                   \n\
                                                                         \n\
	list_mn      : list on ACNUC entry names                         \n\
	CDS_FAM_file : format: FAM <tab> CDS <cr>                        \n\
	outfile      : name of the output file                           \n\
                                                                         \n\n" 
);


}
/********************** end print_mess ********************************/


/*************************** main ***********************************/
main(int argc,char *argv[])
{  
FILE *mn_file;
FILE *outfile;
FILE *cds_file;

if(argc != 4) {
	print_mess();
	exit(1);
   	}

if ((mn_file = fopen(argv[1], "r")) == NULL) {
	fprintf(stderr, "Unable to open file %s\n", argv[1]);
	exit(1);
	}

if ((cds_file = fopen(argv[2], "r")) == NULL) {
	fprintf(stderr, "Unable to open file %s\n", argv[2]);
	exit(1);
	}


if ((outfile = fopen(argv[3], "w")) == NULL) {
	fprintf(stderr, "Unable to create file %s\n", argv[3]);
	exit(1);
	}

acnucopen();
nb_warning = 0;

load_genefam(cds_file);
fclose(cds_file);

do_traitement(mn_file, outfile);

free((char *) cdsfam);

fclose(mn_file);
fclose(outfile);

printf("Normal program end!\n");
exit(0);
}
/********************** end main ********************************/



void load_genefam(FILE *fic)
{
char string[101];
int ii = 0;
char fam[WIDTH_MAX];
char cds[WIDTH_MAX];
char buf[WIDTH_MAX];
int nsub;
off_t faddr;
int div;

nb_seq_tot = 0;

fprintf(stderr, "Loading family / cds file ...   ");
fflush(stderr);


while(fgets(string, 100, fic) != NULL) {
	nb_seq_tot++;
	}
rewind(fic);
fprintf(stderr, "%d sequences\n", nb_seq_tot);
fflush(stderr);

cdsfam = (struct SEQ *)check_alloc(nb_seq_tot,sizeof(struct SEQ));


while(fgets(string, 100, fic) != NULL) {
	if(sscanf(string, "%s %s", fam, cds) != 2) {
		fprintf(stderr, "EXIT: format error (1) \n");
		exit(1);
		}
	if((int) strlen(cds) > WIDTH_MAX) {
		fprintf(stderr, "increase WIDTH_MAX and recompile (1) \n");
		exit(1);
		}

	if((int) strlen(fam) > WIDTH_MAX) {
		fprintf(stderr, "increase WIDTH_MAX and recompile (2) \n");
		exit(1);
		}

		if(! (nsub = isenum(cds))) {
		fprintf(stderr, "EXIT: %s is not in the database !\n",cds);
		exit(1);
		}
		readsub(nsub);
		seq_to_annots64(nsub, &faddr, &div);
		cdsfam[ii].adresse=faddr;
		strcpy(cdsfam[ii].fam, fam);
		strcpy(cdsfam[ii].cds, cds);
		/*cdsfam[ii].adresse=abs(faddr);*/
		Nieme_mot(cds,1,".",buf);
		strcpy(cdsfam[ii].mere, buf);
		if(! (nsub = isenum(buf))) {
		fprintf(stderr, "EXIT: %s is not in the database !\n",cds);
		exit(1);
		}
		readsub(nsub);
		seq_to_annots64(nsub, &faddr, &div);
		cdsfam[ii].madresse=cdsfam[ii].adresse-faddr;
		printf("[%s\t%s\t%lld\t%s\t%lld]\n",cdsfam[ii].fam,cdsfam[ii].cds,cdsfam[ii].adresse,cdsfam[ii].mere,cdsfam[ii].madresse);
		ii++;

	}

/* sort seq. by alphabetical order */
qsort(	cdsfam, 
	nb_seq_tot, 
	sizeof(cdsfam[0]), 
	sort_cdsfam);

printf("Sorted cds:\n");
for (ii=0;ii<nb_seq_tot;ii++){
	printf("[%s\t%s\t%lld\t%lld] ",cdsfam[ii].fam,cdsfam[ii].cds,cdsfam[ii].adresse,cdsfam[ii].madresse);
	if(! (nsub = isenum(cdsfam[ii].cds))) {
	fprintf(stderr, "EXIT: %s is not in the database !\n",cds);
	exit(1);
	}
	readsub(nsub);
	seq_to_annots64(nsub, &faddr, &div);
	printf("\t(%lld) \n",faddr);
	}

}
/********************** end load_cdsfam ********************************/



/********************** sort_cdsfam ********************************/
/* sort seq. by alphabetical order */

int sort_cdsfam(const void *a ,const void *b)
{
struct SEQ *aa, *bb;

aa = (struct SEQ *) a;
bb = (struct SEQ *) b;

/*return (strcmp(aa->cds, bb->cds)); 
return(aa->adresse - bb->adresse);*/
if (aa->adresse > bb->adresse)
		return(1);
if (aa->adresse < bb->adresse)
		return(-1);
if (aa->adresse == bb->adresse)
		return(0);			


}
/********************** end sort_cdsfam ********************************/


void do_traitement(FILE *mn_file, FILE *outfile)
{
char string[101];
char mn[WIDTH_MAX];
int ii = 0;
int nbcds = 0, nb, nbfam = 0;


while(fgets(string, 100, mn_file) != NULL) {
	if(sscanf(string, "%s", mn) != 1) continue;
	nbcds += add_annot(mn, outfile, &nb);
	nbfam += nb;
	ii++;
	}


printf("%d extracted entries.\n", ii);
printf("%d CDS described.\n", nbcds);
printf("%d family written.\n", nbfam);
printf("%d WARNING(S).\n", nb_warning);


}
/********************** end do_traitement ********************************/




int add_annot(char *mn, FILE *outfile, int *nbfam)
{
char fam[WIDTH_MAX];
char cds[WIDTH_MAX];
int test = 0, ncds = 0;
int nsub;
int enr_info;		/* numero d'enr. dans info */
int num_inf;		/* = plinf de subseq pour retrouver les info */
int ii;
off_t faddr,*pfaddr,adr_cds,cur_adr_cds;
int div;
/*off_t nb_bits;*/
int flag_fam,flag_tag;
char ensprot[WIDTH_MAX];
char lamere[WIDTH_MAX],lafille[WIDTH_MAX];


if(! (nsub = isenum(mn))) {
	fprintf(stderr, "EXIT: %s is not in the database !\n", mn);
	exit(1);
	}

readsub(nsub);
num_inf = psub->plinf;

if (psub->pext > 0) {
	fprintf(stderr, "EXIT: sequence %s is a subseq - unexpected!\n", mn);
	exit(1);
	}
 
seq_to_annots64(nsub, &faddr, &div);
read_annots64(faddr, div);
adr_cds=faddr;
/*nb_bits=0;*/
cur_adr_cds=0;

printf("\nCommence  l'adresse %lld \n",adr_cds);
flag_fam=-1;
flag_tag=0;

for(;;) {

	if((strncmp(pinfo->line, "FT   CDS             ", 21) == 0) ||(strncmp(pinfo->line, "FT   mat_peptide     ", 21)) == 0){
	ncds++;
	/*cur_adr_cds=adr_cds+nb_bits;*/
	/*printf("\nNouveau CDS a l'adresse %lld (abs=%lld)\n",cur_adr_cds,abs(cur_adr_cds));*/
	/*printf("\nNouveau CDS a l'adresse %lld ( soit + %lld)\n",cur_adr_cds,cur_adr_cds-adr_cds);
	printf("\nDEBUG: read_annot donne %lld \n",faddr);*/
	cur_adr_cds=faddr;
	printf("\nNouveau CDS a l'adresse %lld ( soit + %lld)\n",cur_adr_cds,cur_adr_cds-adr_cds);
	flag_fam=-1;
	flag_fam=check_fam(cur_adr_cds, outfile);
	/*flag_fam=check_fam(abs(cur_adr_cds), outfile);*/
	printf("resultat de la recherche=%d\n",flag_fam);
	if (flag_fam >=0) {
		printf("Nom du CDS a cette adresse = [%s]\n",cdsfam[flag_fam].cds);
		strcpy(lamere,mn);
		Nieme_mot(cdsfam[flag_fam].cds,1,".",lafille);
		printf("Nom de la mere courante    =             [%s]\n",lamere);
		printf("Nom de la mere de la fille =             [%s]\n",lafille);
		while (strcmp(lafille,lamere) !=0){
			printf("%s n'appartient pas a %s!\nSaute a l'adresse precedente...\n",lafille,lamere);
			flag_fam--;
			Nieme_mot(cdsfam[flag_fam].cds,1,".",lafille);
			printf("Nom de la mere de la fille =             [%s]\n",lafille);
			if (cdsfam[flag_fam].adresse != cur_adr_cds){
				printf("adresse differente!Recherche dans lautre sens\n");
				
	
				
				
				
				break;
					
				}
			}
			while (strcmp(lafille,lamere) !=0){
			printf("%s n'appartient pas a %s!\nSaute a l'adresse suivante...\n",lafille,lamere);
			flag_fam++;
			Nieme_mot(cdsfam[flag_fam].cds,1,".",lafille);
			printf("Nom de la mere de la fille =             [%s]\n",lafille);
			if (cdsfam[flag_fam].adresse != cur_adr_cds){
				printf("adresse differente!\n");
				
				
				
				
				flag_fam=-1;
				break;
					
				}
			}		

		
		if (flag_fam >=0) {
		printf("le CDS %s est selectionne\n",cdsfam[flag_fam].cds);
		flag_tag=1;
		test++;
		} else {
		printf("WARNING! Pas de CDS associé!\n");
		} 
		 
		}
		
						/*
		strcpy(lafille,cdsfam[flag_fam].cds);
		printf("Nom de la mere =             [%s]\n",mn);
		printf("%s n'appartient pas a %s!\nSaute a l'adresse suivante...\n",lafille,lamere);
		fflush(stdout);*/
	/*}
		/*if (strncmp(lafille,lamere,strlen(lamere)) !=0){
		/*while(strncmp(mn,cdsfam[flag_fam].cds,strlen(mn)) !=0) {*/
		/*if (strcmp(lafille,lamere) !=0) {
		printf("%s n'appartient pas a %s!\nSaute a l'adresse suivante...\n",lafille,lamere);
		flag_fam++;	*//*
		exit(1);
		}* /*
		printf("le CDS %s est selectionne\n",cdsfam[flag_fam].cds);
		}*/
		
	}
/*	if(strncmp(pinfo->line, "FT                   /gene=\"", 28) == 0) {
		if (flag_fam >=0) {
	     	fprintf(outfile,  "FT                   /gene_family=\"%s\"\n", cdsfam[flag_fam].fam);
		test++;
		printf("Ecrit la famille no %d : %s (%s)\n",test,cdsfam[flag_fam].fam,cdsfam[flag_fam].cds );
		strcpy(ensprot,cdsfam[flag_fam].cds);
		for (ii=0;ii<strlen(ensprot);ii++) {
			if (ensprot[ii] == '\0') 
					break;
			if (ensprot[ii] == '.') 
				ensprot[ii] = '_';
					
		}
		fprintf(outfile,  "FT                   /db_xref=\"HOMOLENSPROT:%s\"\n", ensprot);
		flag_fam=-1;
	}			
	}
*/		
	fprintf(outfile, "%s\n", pinfo->line);
	
/* attention ne pas ecrire aund on est encore dans le CDS... on choisit le shlash comem signe que c'est ok*/
	
	if (flag_tag) {
		if (strncmp(pinfo->line, "FT                   /", 22) == 0) {
			fprintf(outfile,  "FT                   /cds_name=\"%s\"\n", cdsfam[flag_fam].cds);
			fprintf(outfile,  "FT                   /gene_family=\"%s\"\n", cdsfam[flag_fam].fam);
			flag_tag=0;
			}
		
	}
	/*nb_bits = nb_bits + (off_t) strlen(pinfo->line)+1;*/

	
	if(strncmp(pinfo->line, "//", 2) == 0) break;	
	/*next_annots64(NULL);*/
	next_annots64(&faddr);		
	}


	printf("NOTE: %s  : %d CDS described, %d family written \n", 
		mn, ncds, test);

*nbfam = test;
return(ncds);
}
/********************** end add_annot ********************************/


int check_fam(off_t cur_cds, FILE *outfile)
{

int ii, deb2, fin2, pos, diff;
int test  = 0;

/*printf("Recherche cds %d\n",cur_cds);*/
/* recherche dichotomique 

a modifier si c'est negatif ?*/

deb2 = 0;
fin2 = nb_seq_tot -1;

while(deb2 <= fin2) {
	pos = (deb2 + fin2)/2;
	/*diff = (cur_cds -cdsfam[pos].adresse);*/
	if (cur_cds > cdsfam[pos].adresse)
			diff=1;
	if (cur_cds < cdsfam[pos].adresse)
			diff=-1;
	if (cur_cds == cdsfam[pos].adresse)
			diff=0;	
	/*printf("search %lld %lld %d\n",cur_cds,cdsfam[pos].adresse,diff);*/	
	if(diff<0)  fin2 = pos -1;
	else if(diff>0) deb2 = pos +1;
	else {
		test = 1;
		break;
		}
	}


if(test == 0) {
	return(-1);
	}


return(pos);

}
/********************** end print_fam ********************************/


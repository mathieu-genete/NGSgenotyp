#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <libgen.h>
#include <sys/stat.h>

#define TAILLE_MAX 1000
#define VERSION_NBR "version 1.02"

/*
	compilation: gcc -Wall -Wextra -Werror --std=gnu99 -o splitFastq splitFastq.c
*/

char * rindex(const char *s, int c);

void show_help()
{

	printf("usage : splitFastq [-c cut position] [-o output dir]... Fastq files...\n%s\n\n",VERSION_NBR);
	printf("Split one or more fastq file(s) in 2 files (_R1_RSplit.fastq and _R2_RSplit.fastq)\n\n");
	printf("\t-c\tcut position: if this option is not set the read is half cut\n\n");
	printf("\t-o\toutput directory for splited fastq\n\n");
	printf("\t-h\tshow this help\n\n");
}


int main(int argc, char ** argv)
{
	int manualCut = 0;
	char * outdir = NULL;
	char **fastqFiles;
	int i,c, sizeStr, LenOutName;
	extern char * optarg;
	int chr,id;
	int cutindex = 0;
	char *fname, *tmp;
	char chaine[TAILLE_MAX] = "";
	
	struct stat st;
	
	id = 1;
	
	if (argc==1)
	{
		show_help(argv[0]);
		return 1;
	}
	
	while ((chr = getopt(argc , argv, "hc:o:")) != -1)
	{
		 switch (chr)
		{
			case 'c':
				id+=2;
				cutindex = atoi(optarg);
				manualCut = 1;
				printf("Cutindex = %d\n",cutindex);
				break;
			case 'o':
				id+=2;
				outdir = optarg;
				printf("outdir: %s\n",outdir);
				break;
			case 'h':
				show_help();
				return 1;
			case '?':
				return 1;
		}
	}
	
	if(outdir == NULL || stat(outdir,&st) == 0)
	{
		if (argc>id)
		{
			int fi;
			FILE* fichier = NULL;
			FILE* foutput1 = NULL;
			FILE* foutput2 = NULL;
			FILE* outf = NULL;
		
			fastqFiles = &argv[id];
		
			for(fi=0;fi<(argc-id);fi++)
			{
				printf("-- split: %s --\n",fastqFiles[fi]);
			
				char filename[strlen(fastqFiles[fi])];
			
				strcpy(filename,fastqFiles[fi]);
			
				fichier = fopen(filename, "r+");
			
				fname = basename(filename);
				tmp = rindex(fname,'.');
				tmp[0] = '\0';
			
				LenOutName = strlen(fname)+30;
			
				if (outdir != NULL)
				{
					LenOutName += strlen(outdir);
				}
			
				char out1_name[LenOutName];
				char out2_name[LenOutName];
			
				if (outdir != NULL)
				{
					strcpy(out1_name,outdir);
					strcpy(out2_name,outdir);
					if (outdir[strlen(outdir)-1] != '/')
					{
						strcat(out1_name,"/");
						strcat(out2_name,"/");
					}
					strcat(out1_name,fname);
					strcat(out2_name,fname);
			
				}else
				{
					strcpy(out1_name,fname);
					strcpy(out2_name,fname);
				}
			
				foutput1 = fopen(strcat(out1_name,"_R1_RSplit.fastq"),"w");
				foutput2 = fopen(strcat(out2_name,"_R2_RSplit.fastq"),"w");

				if (fichier != NULL)
				{
					i = 0;
					while (fgets(chaine, TAILLE_MAX, fichier) != NULL)
					{
						c = 0;
						if (manualCut==0)
						{
							sizeStr = strlen(chaine);
							cutindex = sizeStr/2;
						}
					
						while (c<2)
						{	if (c==0)
							{
								outf = foutput2;
							}else
							{
								outf = foutput1;
							}
							switch(i%4)
							{
								case 0:
								  fputs(chaine, outf);
								  break;
								case 1:
								  if (c==1)
								    {
								      chaine[cutindex]='\n';
								      chaine[(cutindex)+1]='\0';
								      fputs(chaine, outf);
								    }else
								    {
								      fputs(&chaine[cutindex], outf);
								    }
								  break;
								case 2:
								  fputs(chaine, outf);
								  break;
								case 3:
								  if (c==1)
								    {
								      chaine[cutindex]='\n';
								      chaine[(cutindex)+1]='\0';
								      fputs(chaine, outf);
								    }else
								    {
								      fputs(&chaine[cutindex], outf);
								    }
								  break;
							}
							c++;
						}
				
						i++;
					}

					fclose(fichier); // On ferme le fichier qui a été ouvert
					printf("%d reads traited for %s\n\n",i/4,fastqFiles[fi]);
				}
			
				fclose(foutput1);
				fclose(foutput2);
			}

		}
	}else
	{
		printf("%s not present\n",outdir);
	}

}

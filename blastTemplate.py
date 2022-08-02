from Bio.Blast.Applications import NcbiblastpCommandline


template_type = 'homo'
froot = 'avastin_5-10mer_0.6_2'
subject = f'templates/{template_type}_templates.fasta'

command = NcbiblastpCommandline(query=f'{froot}/{froot}_modified_sorted.fasta',
                                subject=subject,
                                outfmt=6,
                                out=f'{froot}/{froot}_blast{template_type}Template.m8',
                                # out=f'{froot}/{froot}_blast{template_type}Template.html',
                                # html=True,
                                )
command()



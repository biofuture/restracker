#import Pkg
#Pkg.add("ArgParse")
#Pkg.add("DataStructures")
using DataStructures
using ArgParse
using Dates
import Dates

###
#Developted by Dr. Xiao-Tao Jiang Copyright Hold, March 2019 to March 2021
#Email: biofuture.jiang@gmail.com
#Prerequest software installed
#metaSPAdes
#prodigal
#diamond
###

PWD = pwd()
(SARGDB, SARGSTRUCTURE, MARKER, MPA) = ("$PWD/db/SARG_20170328_12307.dmnd", "$PWD/db/structure_20170328_12037.LIST", "$PWD/db/markers.fasta", "$PWD/db/mpa_v20_m200_marker_info.txt")

function parse_commandline()
	s = ArgParseSettings()

	@add_arg_table s begin
	"--opt1"
		help = "an option with an argument"
	"--opt2"
		help = "another option with an argument"
		arg_type = Int
		default = 0
	"--threads", "-p"
		help = "number of threads used"
		arg_type = Int
		default = 1
	"--memory", "-m"
		help = "memory allocated Gbs"
		arg_type = Int
		default = 1
	"--z"
		help = "an option without argument, i.e. a flag"
		action = :store_true
	"fq1"
		help = "fq 1 file for pair-end reads"
		required = true
	"fq2"
		help = "fq 2 file for pair-end reads"
		required = true
	"--output", "-o"
		help = "Output directory store all results"
		required = true
	end
	return parse_args(s)
end

function main()
	parsed_args = parse_commandline()
	println("Parsed args: $ARGS")
	###
	#for (arg,val) in parsed_args
	#	println(parsed_args[arg])
	#end
	###
	#Perform assemble with metaSPAdes 
	#metaspades.py -1 /srv/scratch/z3524677/Disease_metagenomics/RESISTOME/Host-Tracking/test_1.fq	-2	/srv/scratch/z3524677/Disease_metagenomics/RESISTOME/Host-Tracking/test_2.fq	-t	8	-m	25	-o	/srv/scratch/z3524677/Disease_metagenomics/RESISTOME/Host-Tracking/testout
	#generate the command for metaSPAdes
	#run(`export PATH=\$PATH:$PWD/bin`);

	output = parsed_args["output"]
	fq1 = parsed_args["fq1"]
	fq2 = parsed_args["fq2"]
	threads = parsed_args["threads"]
	memory = parsed_args["memory"]
	assemble_dir = "$output/assemble"
	#run(`perl -e 'unless(-d $output){ mkdir $output }'`);
	if isdir("$output")
	else
		mkdir("$output")
	end 
	if isdir("$assemble_dir")
	else
		mkdir("$assemble_dir")
	end 
	logfile = "$output/log.txt"
	f = open(logfile, "w")	
	time_annotation(f, "Resistome Host Tracking Starting")
	assemble_stdout = "$output/assemble/stdout.log"
	#assemble_command = ["metaspades.py  -1 $fq1 -2 $fq2 -t $threads -m $memory -o $assemble_dir"]
	assemble_command = `metaspades.py  -1 $fq1 -2 $fq2 -t $threads -m $memory -o $assemble_dir`
	#replace(assemble_command, r"\'" => "")
	println("The type info is $assemble_command")
	time_annotation(f, "Assemble pair-end fastq with metaSPAdes\n$assemble_command\t")
	run(pipeline(assemble_command, stdout = "$assemble_stdout"));
	time_annotation(f, "Assemble finished")
	
	#gene prediction with prodigal
	if isdir("$output/gene_prediction")
	else
		mkdir("$output/gene_prediction")
	end 
	
	gene_nuc = "$output/gene_prediction/predict.fna"
	gene_pro = "$output/gene_prediction/predict.faa"
	gene_gbk = "$output/gene_prediction/predict.gbk"	
	contigfa = "$output/assemble/contigs.fasta"
	filtercontig = "$output/filter_length_contig.fasta"
	filter_by_length(contigfa, filtercontig, 100)
	#prodigal -a test.faa -d test.fna -p meta -i outputtest/assemble/contigs.fasta -o test.out -q
	geneprediction_command = `$PWD/bin/prodigal  -a $gene_pro -d $gene_nuc -p meta -i $filtercontig -q -o $gene_gbk`
	time_annotation(f, "Starting gene prediction with Prodigal\n$geneprediction_command\t")
	run(geneprediction_command);
	time_annotation(f, "Finish gene prediction with Prodigal")

	#alignment genes with SARG2.0 database to identify ARGs, we do this alignment with diamond 
	time_annotation(f, "Alignment gene with SARG2.0 database")
	matched_args = "$output/matched_args.fna"
	sarg_alignment = "$output/sarg_alignment.tab.txt"
	argalignment_command = `$PWD/bin/diamond  blastx --db $SARGDB -f 6 --query $gene_nuc --evalue 1e-5 --id 20 --al $matched_args -p $threads -o $sarg_alignment --max-target-seqs 1 --quiet`
	time_annotation(f, "ARGs are identified\n$argalignment_command\t")
	#identify contigs carry ARGs 
	run(argalignment_command);			

	arg_carrier = "$output/arg_carrycontigs.fasta"
	args_neighbourgenes = "$output/arg_neighbourgenes.fasta"
	contig_args_othergenes = "$output/contig_args_othergenes.table.txt"
	argcarry_contigs_genes(filtercontig, sarg_alignment, arg_carrier, gene_nuc, args_neighbourgenes, contig_args_othergenes)	
	
	#alignment with Markers 
	
	#markeralignment_command = `usearch -usearch_global $args_neighbourgenes -db $MARKER -id 0.95 -blast6out $blastout -strand both -maxaccepts 1  -maxrejects 256 -query_cov 0.90`
	#run(markeralignment_command);
	##using minimap2 to maping genes to marker gene sets 
	marker_alignment = "$output/marker_alignment.paf"
	markeralignment_command = `$PWD/bin/minimap2  -c $MARKER  $args_neighbourgenes`
	time_annotation(f, "Algnment other genes with Clade Specific Marker Gene Set\n$markeralignment_command\t")
	run(pipeline(markeralignment_command, stdout="$marker_alignment"));
	
	#taxonomy annotation of contigs gene with StrainPhlan Marker Genes 
	type_strain = "$output/host_strian_arg_type.tab.txt"
	type_spe = "$output/host_species_arg_type.tab.txt"
	type_genus = "$output/host_genus_arg_type.tab.txt"
	subtype_spe = "$output/host_species_arg_subtype.tab.txt"
	subtype_str = "$output/host_strian_arg_subtype.tab.txt"
	subtype_genus = "$output/host_genus_arg_subtypetab.txt"
	detail_contiginfo = "$output/contig_details_arg_taxa.txt"
	gene_subtype_type_carrying_host(markeralignment, MPA, SARGSTRUCTURE, contig_args_othergenes, type_strain, type_spe, type_genus, subtype_spe, subtype_str, subtype_genus, detail_contiginfo)
	#summarize taxonomic profile and the  ARGs subtype and types that carrying 
	
	close(f)
end

##This function process the above results and generate summary information on host carrying resistome 
#

function gene_subtype_type_carrying_host(argneighbours, strianphlan_mpa, sargstructure, contig_arg_neight_map, gene_host, subtype_host, type_host, host_gene, host_subtype, host_type, detail_contiginfo)

         ##contigs taxonomy assignment   contigs => genes => taxonomy  If without annotation, set is as lacking markers as zero 
         contiggenedict = DefaultDict(Dict)   ##using this grammar to define dict's dict 
         matchedmarker = Dict{String, String}() #store all the matched genes in the arg carrying contigs 
         matchedgene = Dict{String, String}()  #store all the  matched markers within strianphlan mpa database 

         contigcov = Dict{String, Float64}()  #store the coverage of each contig 
         contigannotation = Dict{String, String}()   #store the taxonomy of each contig 
         taxacontig = DefaultDict(Dict) ##store all the contigs belong to one taxonomy 

         open(argneighbours) do f
                 line = 1
                 while !eof(f)
                         x = readline(f)
                         items = split(x, "\t")
                         names = split(items[1], "_")
                         contigname = join(names[1:(end-1)], "_")
                         contigcov[contigname] = parse(Float64, names[6])
                         contiggenedict[contigname][items[1]] = items[6]
                         matchedmarker[items[6]] = items[1]
                         matchedgene[items[1]] = items[6]
                         line += 1
                 end
         end

        ##store all the names of markers 
         markerstaxa = Dict{String, String}()  ##store the last level of taxonomy of the mateched markers
         markersphylo = Dict{String, String}() ##store the whole taxonomy lineage of the matched markers 
         open(strianphlan_mpa) do f
                 line = 1
                 while !eof(f)
                         x = readline(f)
                         items = split(x, "\t")

                         names = match(r"\'clade\': (\S+)\, ", items[2]).match
                         phylo = match(r"\'taxon\': (\S+)\}", items[2]).match
                         names = (split(names, ": "))[2]
                         phylo = (split(phylo, ": "))[2]
                         names = replace(names, r", |'" => "")
                         phylo = replace(phylo, r"}|'" => "")
                         if  haskey(matchedmarker, items[1]) == true
                                 markerstaxa[items[1]] =  names
                                 markersphylo[items[1]] = phylo
                                 #println("$names\t$phylo")
                         end
                         line += 1
                 end
         end
         ##store subtype and type information of args  
         genesubtype = Dict{String, String}()  ##store all the subtype information of each arg 
         genetype = Dict{String, String}()     ##store the type information of each arg 

         open(sargstructure) do f
                 line = 1
                 while !eof(f)
                         x = readline(f)
                         items = split(x, "\t")
                         items[2] = replace(items[2], r"\[|\]|\'" => "")
                         names = split(items[2], ", ")
                         for na in names
                                genesubtype[na] = items[1]
                                tyname = (split(items[1], "__"))[1]
                                genetype[na] = tyname
                                #println("$na")                         
                         end
                         line += 1
                 end
         end

      ##give each matched contig an taxonomy annotation 
        for (key, val) in  contiggenedict

                #rown = length(keys(contiggenedict[key])) #number of genes carried by this contig 
                #println(keys(contiggenedict[key]))
                arraytaxa = fill("undef", (100,8)) ##store all the taxa rank into this two dimension array; allow a maximum of 100 genes                
                #println(arraytaxa[1, 2])

                i = 1 ##this is index for taxonomy 
                for (k2, v2) in contiggenedict[key]
                        ##store all the markers into dict 
                        if haskey(markerstaxa, v2) == true && haskey(markersphylo, v2) == true
                                #println("$key\t", markerstaxa[v2], "\t",markersphylo[v2])
                                #consunsus_taxa{key}{v2} = markersphylo[v2]
                                separray = split(markersphylo[v2], "|")
                                #println(typeof(separray), "\t", markersphylo[v2], "\t", 1:length(separray),"\t", rown)
                                for j in 1:length(separray)
                                        arraytaxa[i, j] = separray[j]
                                end
                                i += 1
                        end
                        ##consensus taxonomy based on the clade-specific markers 
                end

                ##using lca to generate the consensus taxa 
                #processing the arrartaxa 
                rows = i-1
                #println(rows)
                taxarray = fill("undef", 8)
                for j in 1:8
                        #println(i)
                        dicarray = Dict{String, Int64}()
                        for i in 1:rows
                                #println(arraytaxa[i, j]) 
                                if haskey(dicarray, arraytaxa[i, j])
                                        dicarray[arraytaxa[i, j]] += 1
                                else
                                        dicarray[arraytaxa[i, j]] = 1
                                end
                        end

                        ##sort dicarray by value 
                        asarray = collect(dicarray)
                        assortedvector = sort(asarray, by=x->last(x))
                        consenslevel = "mark"  # mark to indicate successful 
                        for k in 1:length(assortedvector)

                                if  (assortedvector[k][2] / rows) >= 0.667  && assortedvector[k][1] != "undef" ##if over 2/3 of the taxonomy are consistent we give it the taxonomy assignment 
                                        consenslevel = assortedvector[k][1]
                                end
                        end

                        if consenslevel != "mark"
                                taxarray[j] = consenslevel
                        else
                                #taxarray[j] = "missing"
                                @goto escape_label
                        end

                end
                @label escape_label
		#B
                #subarray = deleteat!(taxarray, taxarrau == "undef")
                subarray = filter!(x->x≠"undef",taxarray)
                taxaannotate = join(subarray, "|")  ##consensus taxonomy of these genes annotated with strianphlan marker genes 
                #println(taxaannotate)
                contigannotation[key] = taxaannotate
                taxacontig[taxaannotate][key] = 1
        end
      #
        fcontiginfo = open(detail_contiginfo, "w")
        contigarg = DefaultDict(Dict)
        open(contig_arg_neight_map) do f
                line = 1
                ##process each arg carrying contig to get their args and taxonomy information 
                while !eof(f)
                        x = readline(f)
                        items = split(x, "\t")
                        len = length(items)
                        indexsep = findfirst(x -> x=="OTHERGENE", items)
                        #println(items[3], "is here")   
                        #starting from 4 
                        for i in 4:(indexsep-1)
                                geneargs = split(items[i], ":")
                                if haskey(genesubtype, geneargs[2]) == true
                                        #println(items[1], "\t", geneargs[1], "\t",geneargs[2], "\t",genesubtype[geneargs[2]], "\t",genetype[geneargs[2]])
                                        if haskey(contigarg[items[1]], geneargs[2]) == true
                                                contigarg[items[1]][geneargs[2]] += 1  #copies of this gene on the contig items[1]
                                        else
                                                contigarg[items[1]][geneargs[2]] = 1
                                        end
                                end
                        end

                        if len >  indexsep 
                                for i in (indexsep+1):len
                                        if haskey(matchedgene, items[i]) == true
                                                gid = matchedgene[items[i]]
                                                if haskey(markerstaxa, gid) == true && haskey(markersphylo, gid) == true
							infocontig  = items[1] * items[i] * gid * markerstaxa[gid] * markersphylo[gid]
							write(fcontiginfo, "$infocontig\n")
                                                end
                                        end
                                end
                        end 

                        line += 1
                end
        end

        ##summerize information 

        taxagene = DefaultDict(Dict)   ##including all contigs annotatoin 
        taxasubtype = DefaultDict(Dict)
        taxatype = DefaultDict(Dict)

        alltypes = Dict{String, Int}()
        allgenes = Dict{String, Int}()
        allsubtypes = Dict{String, Int}()


        striangene = DefaultDict(Dict)  ##only including arg carrying contigs that annotate to strain level 
        striansubtype = DefaultDict(Dict)
        striantype = DefaultDict(Dict)

        speciesgene = DefaultDict(Dict)  ##only including arg carrying contigs that annotate to species level 
        speciessubtype = DefaultDict(Dict)
        speciestype = DefaultDict(Dict)

        genusgene = DefaultDict(Dict)  ##only including arg carrying contigs that annotate to genus level (include those to species level)
        genussubtype = DefaultDict(Dict)
        genustype = DefaultDict(Dict)

        ##For the length of args 
        for (taxa, v1) in taxacontig
                for (contig, v2) in taxacontig[taxa]   ##parse every contigs annotate to this taxa 
                        if haskey(contigarg, contig) == true    ##if it is arg carrying contig
                                for (gene, v) in contigarg[contig]        ##get the arg in that contig and their copies v 
                                        if haskey(contigcov, contig) == true   ##get the coverage of this contig

                                                taxagene[taxa][gene] = v * contigcov[contig]  ##sum up arg gene information  ## v is the copies of this arg gene in this contig, usually it should be 1
                                                allgenes[gene] =1
                                                if haskey(genesubtype, gene) == true   #sumup subtye coverage 
                                                        allsubtypes[genesubtype[gene]] = 1
                                                        if haskey(taxasubtype[taxa], genesubtype[gene]) == true
                                                                taxasubtype[taxa][genesubtype[gene]] += contigcov[contig] * v
                                                        else
                                                                taxasubtype[taxa][genesubtype[gene]] = contigcov[contig] * v
                                                        end
                                                end

                                                if haskey(genetype, gene) == true   #sumup type coverage 
                                                        alltypes[genetype[gene]] = 1
                                                        if haskey(taxatype[taxa], genetype[gene]) == true
                                                                taxatype[taxa][genetype[gene]] += contigcov[contig] * v
                                                        else
                                                                taxatype[taxa][genetype[gene]] = contigcov[contig] * v
                                                        end
                                                end


                                                ##get the species information ready 
                                                arraytax = split(taxa, "|")
                                                if length(arraytax) == 8  ##to strain level remove the strain information, as strain information 

                                                        striangene[taxa][gene]= v * contigcov[contig]

                                                        if haskey(genesubtype, gene) == true   #sumup subtye coverage 
                                                                if haskey(striansubtype[taxa], genesubtype[gene]) == true
                                                                        striansubtype[taxa][genesubtype[gene]] += contigcov[contig] * v
                                                                else
                                                                        striansubtype[taxa][genesubtype[gene]] = contigcov[contig] * v
                                                                end
                                                        end

                                                        if haskey(genetype, gene) == true   #sumup type coverage 

                                                                if haskey(striantype[taxa], genetype[gene]) == true
                                                                        striantype[taxa][genetype[gene]] += contigcov[contig] * v
                                                                else
                                                                        striantype[taxa][genetype[gene]] = contigcov[contig] * v
                                                                end
                                                        end

                                                        getspecies = join(arraytax[1:7], "|")

                                                        ##init these species that have the strian information 
                                                        speciesgene[getspecies][gene]= v * contigcov[contig]

                                                        if haskey(genesubtype, gene) == true   #sumup subtye coverage 
                                                                if haskey(speciessubtype[getspecies], genesubtype[gene]) == true
                                                                        speciessubtype[getspecies][genesubtype[gene]] += contigcov[contig] * v
                                                                else
                                                                        speciessubtype[getspecies][genesubtype[gene]] = contigcov[contig] * v
                                                                end
                                                        end

                                                        if haskey(genetype, gene) == true   #sumup type coverage 

                                                                if haskey(speciestype[getspecies], genetype[gene]) == true
                                                                        speciestype[getspecies][genetype[gene]] += contigcov[contig] * v
                                                                else
                                                                        speciestype[getspecies][genetype[gene]] = contigcov[contig] * v
                                                                end
                                                        end


                                                        getgenus = join(arraytax[1:6], "|")

                                                        ##init these species that have the strian information 
                                                        genusgene[getgenus][gene]= v * contigcov[contig]

                                                        if haskey(genesubtype, gene) == true   #sumup subtye coverage 
                                                                if haskey(genussubtype[getgenus], genesubtype[gene]) == true
                                                                        genussubtype[getgenus][genesubtype[gene]] += contigcov[contig] * v
                                                                else
                                                                        genussubtype[getgenus][genesubtype[gene]] = contigcov[contig] * v
                                                                end
                                                        end

                                                        if haskey(genetype, gene) == true   #sumup type coverage 

                                                                if haskey(genustype[getgenus], genetype[gene]) == true
                                                                        genustype[getgenus][genetype[gene]] += contigcov[contig] * v
                                                                else
                                                                        genustype[getgenus][genetype[gene]] = contigcov[contig] * v
                                                                end
                                                        end



                                                elseif length(arraytax) == 7


                                                        ##init these species that have the strian information 
                                                        if haskey(speciesgene[taxa], gene) == true
                                                                speciesgene[taxa][gene] += v * contigcov[contig]
                                                        else
                                                                speciesgene[taxa][gene] = v * contigcov[contig]
                                                        end
                                                        if haskey(genesubtype, gene) == true   #sumup subtye coverage 
                                                                if haskey(speciessubtype[taxa], genesubtype[gene]) == true
                                                                        speciessubtype[taxa][genesubtype[gene]] += contigcov[contig] * v
                                                                else
                                                                        speciessubtype[taxa][genesubtype[gene]] = contigcov[contig] * v
                                                                end
                                                        end

                                                        if haskey(genetype, gene) == true   #sumup type coverage 

                                                                if haskey(speciestype[taxa], genetype[gene]) == true
                                                                        speciestype[taxa][genetype[gene]] += contigcov[contig] * v
                                                                else
                                                                        speciestype[taxa][genetype[gene]] = contigcov[contig] * v
                                                                end
                                                        end


                                                        getgenus = join(arraytax[1:6], "|")

                                                        ##init these species that have the strian information 
                                                        if haskey(genusgene[getgenus], gene) == true
                                                                genusgene[getgenus][gene] += v * contigcov[contig]
                                                        else
                                                                genusgene[getgenus][gene] = v * contigcov[contig]
                                                        end

                                                        if haskey(genesubtype, gene) == true   #sumup subtye coverage in the genus which can be assinged to a species 
                                                                if haskey(genussubtype[getgenus], genesubtype[gene]) == true
                                                                        genussubtype[getgenus][genesubtype[gene]] += contigcov[contig] * v
                                                                else
                                                                        genussubtype[getgenus][genesubtype[gene]] = contigcov[contig] * v
                                                                end
                                                        end

                                                        if haskey(genetype, gene) == true   #sumup type coverage 

                                                                if haskey(genustype[getgenus], genetype[gene]) == true
                                                                        genustype[getgenus][genetype[gene]] += contigcov[contig] * v
                                                                else
                                                                        genustype[getgenus][genetype[gene]] = contigcov[contig] * v
                                                                end
                                                        end

                                                elseif length(arraytax) == 6

                                                ##get the genus information ready 
                                                        ##init these species that have the strian information 
                                                        if haskey(genusgene[taxa], gene) == true
                                                                genusgene[taxa][gene] += v * contigcov[contig]
                                                        else
                                                                genusgene[taxa][gene] = v * contigcov[contig]
                                                        end
                                                        if haskey(genesubtype, gene) == true   #sumup subtye coverage 
                                                                if haskey(genussubtype[taxa], genesubtype[gene]) == true
                                                                        genussubtype[taxa][genesubtype[gene]] += contigcov[contig] * v
                                                                else
                                                                        genussubtype[taxa][genesubtype[gene]] = contigcov[contig] * v
                                                                end
                                                        end

                                                        if haskey(genetype, gene) == true   #sumup type coverage 

                                                                if haskey(genustype[taxa], genetype[gene]) == true
                                                                        genustype[taxa][genetype[gene]] += contigcov[contig] * v
                                                                else
                                                                        genustype[taxa][genetype[gene]] = contigcov[contig] * v
                                                                end
                                                        end
                                                else
                                                        ##not count this part of contigs, they are with too high level of taxonomy assignment 
                                                end
                                        end
                                end
                        end

                end
        end

        ##output taxatype 
        sortedalltypes = sort(collect(keys(alltypes)))
        sortedallgenes = sort(collect(keys(allgenes)))
        sortedallsubtypes = sort(collect(keys(allsubtypes)))

	#gene_subtype_type_carrying_host(argneighbours, strianphlan_mpa, sargstructure, contig_arg_neight_map, gene_host, subtype_host, type_host, host_gene, host_subtype, host_type, detail_contiginfo)


        fstrian = gene_host
        fspe = subtype_host
        fgenus = type_host
        spesubtype = host_gene
        strsubtype = host_subtype
        gensubtype = host_type 

        write_out(fspe,speciestype,sortedalltypes)
        write_out(fstrian,striantype,sortedalltypes)
        write_out(fgenus,genustype,sortedalltypes)
        write_out(spesubtype, speciessubtype, sortedallsubtypes)
        write_out(strsubtype, striansubtype, sortedallsubtypes)
        write_out(gensubtype, genussubtype, sortedallsubtypes)

end

function write_out(fout, taxaleveldict, sorttypesarray)

        fo = open(fout, "w")
        typesheader = join(sorttypesarray, "\t")
        #println("$typesheader")
        write(fo, "\t$typesheader\n")
        for (tax, v1) in taxaleveldict
                        if(tax == "")
                                tax = "Unclassified"
                        end
                        line = tax
                        for ty in sorttypesarray
                #                       println(line)
                                if haskey(taxaleveldict[tax], ty) == true
                                        line = line * "\t" * string(taxaleveldict[tax][ty])
                                else
                                        line = line * "\t" * string("0")
                                end

                        end
                        write(fo, "$line\n")
        end


end

##host carrying tree
##this  function parse the alignment results and the contigs that carry  args
function argcarry_contigs_genes(inputcontig, sarg_alignment, arg_carriercontig, genein, geneout, contig_args_othergenes)

	contigarg = DefaultDict(Dict)
	contigother = DefaultDict(Dict)
        
	sargdict = Dict{String, String}()
        contigdict = Dict{String, String}()

        open(sarg_alignment) do f
                line = 1
                while !eof(f)
                        x = readline(f)
                        tabarray = split(x, "\t")
                        sargdict[tabarray[1]] = tabarray[2]

                        namearray = split(tabarray[1], "_")
                        pop!(namearray)
                        contign = join(namearray, "_")
                        contigdict[contign] = "1"

			contigarg[contign][tabarray[1]] = tabarray[2]
                        line += 1
                end
        end

        selectfasta_by_name(inputcontig, contigdict, arg_carriercontig)
        ##fetch neighbourgenes

        finput = read(genein, String)
        foutput = open(geneout, "w")
        fileblocks = split(finput, ">")
        len = length(fileblocks)
	
        for i in 2:len
                nameseq = split(fileblocks[i], "\n")
                name = split(nameseq[1], r"\s+")
                namecontig = split(name[1], '_')
                pop!(namecontig)
                namecontig = join(namecontig, '_')

                deleteat!(nameseq, 1)
                deleteat!(nameseq, length(nameseq))
                seq = join(nameseq, "\n")
                #println(name[1])
                if (haskey(contigdict, namecontig) == true)  &&  (haskey(sargdict, name[1]) == false)
                        write(foutput, ">",name[1], "\n", seq,"\n")
        		contigother[namecontig][name[1]] = 1
		end
        end

	##ftab record mapping information of contig carraied args and neighour genes 
	ftab  = open(contig_args_othergenes, "w")		
	for  (cg, value) in contigdict
		argstring = ""
		othergenestring = ""
		if haskey(contigarg,  cg) == true
			for (ke, va) in contigarg[cg]	
				s = join([ke, va], ":")
				argstring = join([argstring, s], "\t")		
			end 
		end 

		if (haskey(contigother, cg) == true)
			for (ke, va) in contigother[cg]	
				othergenestring = join([othergenestring, ke], "\t")		
			end 
			
		end 
		write(ftab, "$cg\tARGS\t$argstring\tOTHERGENE\t$othergenestring\n")
	end

end


function selectfasta_by_name(inputfa, dict, outputfa)
	#dict include all the names of contigs that needed 
	finput = read(inputfa, String)
	foutput = open(outputfa, "w")
	fileblocks = split(finput, ">")
	len = length(fileblocks)

	for i in 2:len
		nameseq = split(fileblocks[i], "\n") 
		name = nameseq[1] 
		deleteat!(nameseq, 1)        
		deleteat!(nameseq, length(nameseq))        
		seq = join(nameseq, "\n")
	
		if haskey(dict, name) == true
			write(foutput, ">$name", "\n", seq,"\n");
		end
	end
end


function filter_by_length(inputfa, outputfa, leng)

 	finput = read(inputfa, String) 
        foutput = open(outputfa, "w")
        #As for one sample the contig file won't be too big 
        fileblocks = split(finput, ">")
        len = length(fileblocks)
	#println("The number of blocks are $len and the length cut of is $leng")
        for i in 2:len
                nameseq = split(fileblocks[i], "\n") 
		name = nameseq[1] 
		deleteat!(nameseq, 1)        
		deleteat!(nameseq, length(nameseq))        
		seq = join(nameseq, "\n")
		#println("$name\n$seq")
		#return
                namearray = split(name, "_" )         
		#t1 = typeof(namearray[4])
		#t2 = typeof(leng)
		num1 = parse(Int64,namearray[4])
		num2 = leng
		#println("The numbers are $t1 and $t2")
		if  num1 >= num2
			#println(namearray[4], "\tis bigger than\t", leng)
                        write(foutput, ">$name", "\n", seq,"\n");
			#println("$name  is not meet the length requirement")
                end                     
        end
        close(foutput)
end

function time_annotation(f, info)
	#write the time and relevent program information to log.txt file
	timenow = Dates.now()
	write(f, "$info: $timenow\n")
end 

main()

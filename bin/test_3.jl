#import Pkg
#Pkg.add("DataStructures")
using DataStructures

function gene_subtype_type_carrying_host(argneighbours, strianphlan_mpa, sargstructure, contig_arg_neight_map, gene_host, subtype_host, type_host, host_gene, host_subtype, host_type)

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
		#subarray = deleteat!(taxarray, taxarrau == "undef")
		subarray = filter!(x->x≠"undef",taxarray)
		taxaannotate = join(subarray, "|")  ##consensus taxonomy of these genes annotated with strianphlan marker genes 
		#println(taxaannotate)
		contigannotation[key] = taxaannotate
		taxacontig[taxaannotate][key] = 1
	end 
	#
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
"""
			if len >  indexsep 
				for i in (indexsep+1):len
					if haskey(matchedgene, items[i]) == true
						gid = matchedgene[items[i]]
						if haskey(markerstaxa, gid) == true && haskey(markersphylo, gid) == true
							println(items[1], "\t",items[i], "\t", gid, "\t",markerstaxa[gid], "\t",markersphylo[gid])
						end
					end
				end
			end 
"""
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
			if haskey(contigarg, contig) == true	##if it is arg carrying contig
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

	typesheader = join(sortedalltypes, "\t")
	println("\t$typesheader")
"""
#	for (taxa, v1) in taxatype
#			if(taxa == "")
#				taxa = "Unclassified"
#			end
#			print(taxa)
#			for ty in sortedalltypes
#				if haskey(taxatype[taxa], ty) == true
#					print("\t",taxatype[taxa][ty])
#				else
#					print("\t0")
#				end
#
#			end
#			println()
#		
#	end 
"""

	for (species, v1) in speciestype
			if(species == "")
				species = "Unclassified"
			end
			print("$species")
			for ty in sortedalltypes
				if haskey(speciestype[species], ty) == true
					print("\t",speciestype[species][ty])
				else
					print("\t0")
				end

			end
			println()
		
	end 
end

gene_subtype_type_carrying_host(ARGS[1], ARGS[2], ARGS[3], ARGS[4], ARGS[5], ARGS[6], ARGS[7], ARGS[8], ARGS[9], ARGS[10])




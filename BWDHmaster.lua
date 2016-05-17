--Start: these values can be changed between 0 and 2
MutationConnectionsChance = 0.25
PerturbChance = 0.90
CrossoverChance = 0.75
LinkMutationChance = 2.0
NodeMutationChance = 0.50
BiasMutationChance = 0.40
StepSize = 0.1
DisableMutationChance = 0.4
EnableMutationChance = 0.2
--end

--This variable can be changed to any number greater than 1.
MaxPopulation = 100

--This variable controls how long the program will wait for mario to move
TimeoutConstant = 20

--This variable is used in the fitness function to remove species that aren't improving
StaleSpecies = 15

--This variable is the save state used to start each run
Filename = "DP1.state"
-- change this to true to output test data to a text file
testDataOn = false
-- modify this variable when running tests with mutations
Description = "Basic Run"
-- set this to the name of your .txt file (will be created in the same folder as the lua script)
DatafileName = "basicRun.txt"
-- if you want to load a file add txt here
-- example: "backup.3.level1testTopFitness.txt" for the 3rd gen backup where DatafileName = "level1testTopFitness.txt"
loadFileName = ""

--Start: don't change these global variables
DeltaDisjoint = 2.0
DeltaWeights = 0.4
DeltaThreshold = 1.0
GlobalName = 0
MaxNeurons = 1000000
BoxRadius = 6
InputSize = (BoxRadius*2+1)*(BoxRadius*2+1)
Inputs = InputSize+1
--End

--setting the buttons available to choose from
ButtonNames = {
"A",
"B",
"X",
"Y",
"Up",
"Down",
"Left",
"Right",
}

Outputs = #ButtonNames --returns the length of the ButtonNames (8)

---------------------------------------------------------------------------------------------------
-- Creating basic populations, species, and chromosomes -------------------------------------------
-- Created: March 4, 2016
-- Last modified: April 6, 2016
-- Primary author: Christopher Bruns

--creates a new population of species
function newPopulation()
  local population = {}
  population.species = {} --contains table of species within the population
  population.generation = 0
  population.innovations = Outputs
  population.currentSpecies = 1
  population.currentChromosome = 1
  population.currentFrame = 0
  population.maxFitness = 0

  return population
end

--creates a new species that is part of a population
function newSpecies()
  local species = {}
  species.topFitness = 0
  species.staleness = 0
  species.chromosomes = {} --contains a table of chromosomes within the species
  species.averageFitness = 0

  return species
end

--creates a new chromosome that is part of a species
function newChromosome()
  local chromosome = {}
  chromosome.genes = {} --has a table of it's genes within the chromosome
  chromosome.fitness = 0
  chromosome.adjustedFitness = 0
  chromosome.network = {} --has information regarding the neural network and links
  chromosome.maxneuron = 0
  chromosome.globalRank = 0
  chromosome.name = GlobalName
  GlobalName = GlobalName + 1
  chromosome.generation = 0
  chromosome.mutationRates = {} --a table containing the rates of different mutations to be taken on the chromosome
  chromosome.mutationRates["node"] = NodeMutationChance
  chromosome.mutationRates["step"] = StepSize
  chromosome.mutationRates["enable"] = EnableMutationChance
  chromosome.mutationRates["disable"] = DisableMutationChance
  chromosome.mutationRates["link"] = LinkMutationChance
  chromosome.mutationRates["bias"] = BiasMutationChance
  chromosome.mutationRates["connections"] = MutationConnectionsChance

  return chromosome
end

--creates a copy of a specific chromosome, copying the genes and the mutationRates
function copyChromosome(chromosome)
  local chromosome2 = newChromosome()
  for g=1,#chromosome.genes do
    table.insert(chromosome2.genes, copyGene(chromosome.genes[g]))
  end
  chromosome2.maxneuron = chromosome.maxneuron
  chromosome2.mutationRates["connections"] = chromosome.mutationRates["connections"]
  chromosome2.mutationRates["link"] = chromosome.mutationRates["link"]
  chromosome2.mutationRates["bias"] = chromosome.mutationRates["bias"]
  chromosome2.mutationRates["node"] = chromosome.mutationRates["node"]
  chromosome2.mutationRates["enable"] = chromosome.mutationRates["enable"]
  chromosome2.mutationRates["disable"] = chromosome.mutationRates["disable"]
  chromosome2.mutationRates["step"] = chromosome.mutationRates["step"]

  return chromosome2
end

--creates a basic chromosome and also calls for a mutation of that chromosome
function basicChromosome()
  local chromosome = newChromosome()
  local innovation = 1

  chromosome.maxneuron = Inputs
  mutateChromosome(chromosome)

  return chromosome
end

---------------------------------------------------------------------------------------------------
-- Creating basic genes ---------------------------------------------------------------------------
-- Created: March 7, 2016
-- Primary author: Christopher Bruns

--creates a new gene, sets all the attributes to initial values and enables the gene
function newGene()
  local gene = {}
  gene.into = 0
  gene.out = 0
  gene.weight = 0.0
  gene.enabled = true
  gene.innovation = 1

  return gene
end

--Copies an existing gene copying all of the values into the copy and returns the copy
function copyGene(gene)
  local gene2 = newGene()
  gene2.into = gene.into
  gene2.out = gene.out
  gene2.weight = gene.weight
  gene2.enabled = gene.enabled
  gene2.innovation = gene.innovation

  return gene2
end


---------------------------------------------------------------------------------------------------
-- Creating neurons, links, and the neural network ------------------------------------------------
-- Created: March 5, 2016
-- Modified: April 6, 2016
-- Primary author: Jack Van Gent


--creats a neuron to be used with the genes within chromosomes
function createNeuron()
  local neuron = {}
  neuron.incomingValues = {}
  neuron.value = 0.0
  return neuron
end

--selecting a random neuron
function randomNeuron(genes, nonInput)
--creating a table of neurons
  local neurons = {}
  --checking if the inputs is empty before adding to the table of neurons
  if not nonInput then
    for i=1,Inputs do
      neurons[i] = true
    end
  end
  --for each item in the output table (length 8) setting it to true (active)
  for o=1,Outputs do
    neurons[MaxNeurons+o] = true
  end
  --[[for each gene in the gene table if it has inputs OR if it has information/neuron link coming into it and that is
  greater than the Inputs then it is updated in the neurons table created above as well. Then do the same for the
  out value as well.]]
  for i=1,#genes do
    if (not nonInput) or genes[i].into > Inputs then
      neurons[genes[i].into] = true
    end
    if (not nonInput) or genes[i].out > Inputs then
     neurons[genes[i].out] = true
    end
  end
  --[[find the total number of neurons that were added to the table from the gene table. Then will use the random
  function to pick out one of them and that is referenced as 'n']]
  local count = 0
  for _,_ in pairs(neurons) do
    count = count + 1
  end
  local n = math.random(1, count)

  for k,v in pairs(neurons) do
    n = n-1
    if n == 0 then
      return k
    end
  end

  return 0
end

--builds the neural network from the chromosome's genes
function createNetwork(chromosome)
  local network = {}
  network.neurons = {}

  --create a neuron for each possible input, should create 170 of these
  for i=1,Inputs do
    network.neurons[i] = createNeuron()
  end

  --create a neuron for each possible output (button press), should create 8 of these
  for j=1,Outputs do
    network.neurons[MaxNeurons+j] = createNeuron()
  end

  table.sort(chromosome.genes, function (a,b)
    return (a.out < b.out)
  end)
  --[[creating neurons based on gene.out and gene.into values to determine if they exist to start
  with and then adding those new neurons]]
  for i=1,#chromosome.genes do
    local gene = chromosome.genes[i]
    if gene.enabled then
      if network.neurons[gene.out] == nil then
        network.neurons[gene.out] = createNeuron()
      end
      local neuron = network.neurons[gene.out]
      table.insert(neuron.incomingValues, gene)
      if network.neurons[gene.into] == nil then
        network.neurons[gene.into] = createNeuron()
      end
    end
  end
  chromosome.network = network
end

--increment the number of innovations within a population
function newInnovation()
	population.innovations = population.innovations + 1
	return population.innovations
end


---------------------------------------------------------------------------------------------------
-- Crossover techniques used during breeding of chromosomes ---------------------------------------
-- Created: March 9, 2016
-- Modified: April 15, 2016
-- Primary author: Christopher Bruns


--first crossover mutation technique
function crossoverMutation(chromosome1, chromosome2)
  --first set chromosome1 to be the one with the best fitness for ease when dealing w/ next sections
  if chromosome2.fitness > chromosome1.fitness then
    tempChrom = chromosome1
    chromosome1 = chromosome2
    chromosome2 = tempChrom
  end

  --creating new child to be made by the two chromosomes crossover
  local childChromosome = newChromosome()

  --creating the innovations table to be used with the genes with the new childChromosome
  local innovations2 = {}

  --Cycle through the length of the table in chromosome2 creating genes for each item in the
  --table and then adding those genes to the innovations2 table.
  for i = 1, #chromosome2.genes do
    local gene = chromosome2.genes[i]
    innovations2[gene.innovation] = gene
  end

  --Now it goes through chromosome1's genes.
  for i = 1, #chromosome1.genes do
    local gene1 = chromosome1.genes[i]
    local gene2 = innovations2[gene1.innovation]
    --Checking to see if gene2 is NOT null and then random 1 or 2, if it returns 1 and it is not
    --null it checks if the gene is enabled, if it is gene2 is placed in the new childChromosome
    --otherwise gene1 is placed within the child
    if gene2 ~= nil and math.random(2) == 1 and gene2.enabled then
      table.insert(childChromosome.genes, copyGene(gene2))
    else
      table.insert(childChromosome.genes, copyGene(gene1))
    end
  end

  --Setting the maxneuron within the child now based on the max of chromosome 1 and 2
  childChromosome.maxneuron = math.max(chromosome1.maxneuron, chromosome2.maxneuron)

  --Takes the type of mutation from the mutationRates table and the rate associated with
  --each type and applies them to the childChromosome that has been created.
  for mutation, rate in pairs(chromosome1.mutationRates) do
    childChromosome.mutationRates[mutation] = rate
  end

  --Return the newly created chromosome from the crossover
  return childChromosome


end

--second crossover mutation technique
function onePointCrossover(chromosome1, chromosome2)
  if chromosome2.fitness > chromosome1.fitness then
    tempChrom = chromosome1
    chromosome1 = chromosome2
    chromosome2 = tempChrom
  end

  --creating new child to be made by the two chromosomes crossover
  local childChromosome = newChromosome()

  --creating the innovations table to be used with the genes with the new childChromosome
  local innovations2 = {}

  for i = 1, #chromosome2.genes do
    local gene = chromosome2.genes[i]
    innovations2[gene.innovation] = gene
  end

  --determining the crossover split point and checking it is a valid split point for both chromosomes
  local split = math.random(1,#chromosome2.genes)
  if split >= #chromosome1.genes then split = (#chromosome1.genes - 1) end

  --[[sections below inserting genes into the child until hitting the split then doing it
  from them other chromosome]]
  for i = 1, split do
    local gene = chromosome1.genes[i]
    if gene ~= nil then table.insert(childChromosome.genes, copyGene(gene)) end
  end

  for i = split, #chromosome2.genes do
    local gene = chromosome2.genes[i]
    if gene ~= nil then table.insert(childChromosome.genes, copyGene(gene)) end
  end

  --setting the maxneuron and the mutationRates for the child.
  childChromosome.maxneuron = math.max(chromosome1.maxneuron, chromosome2.maxneuron)

  for mutation, rate in pairs(chromosome1.mutationRates) do
    childChromosome.mutationRates[mutation] = rate
  end

  --returns only one child due to sorting of the chromosomes according to fitness at the beginning of the function
  return childChromosome
end

--third crossover mutation technique
function twoPointCrossover(chromosome1, chromosome2)
  if chromosome2.fitness > chromosome1.fitness then
    tempChrom = chromosome1
    chromosome1 = chromosome2
    chromosome2 = tempChrom
  end

  --creating new children to be made by the two chromosomes crossover
  local childChromosome1 = newChromosome()
  local childChromosome2 = newChromosome()
  local split1 = math.random(1,#chromosome2.genes)
  local split2 = math.random(1,#chromosome2.genes)

  --setting up two different split points and comparing to make sure works for chromosome1 as well
  if split1 > split2 then
    tempNum = split1
    split1 = split2
    split2 = tempNum
  end
  if split2 > #chromosome1.genes then
    split2 = (#chromosome1.genes - 1)
  end

  --[[based on the two splits checking for nil genes and inserting the genes from the
  parents into the children based on the splits for the two children]]
  for i = 1, split1 do
    local gene = chromosome1.genes[i]
    if gene ~= nil then table.insert(childChromosome1.genes, copyGene(gene)) end
  end
  for i = split1, split2 do
    local gene = chromosome2.genes[i]
    if gene ~= nil then table.insert(childChromosome1.genes, copyGene(gene)) end
  end
  for i = split2, #chromosome1.genes do
    local gene = chromosome1.genes[i]
    if gene ~= nil then table.insert(childChromosome1.genes, copyGene(gene)) end
  end
  for i = 1, split1 do
    local gene = chromosome2.genes[i]
    if gene ~= nil then table.insert(childChromosome2.genes, copyGene(gene)) end
  end
  for i = split1, split2 do
    local gene = chromosome1.genes[i]
    if gene ~= nil then table.insert(childChromosome2.genes, copyGene(gene)) end
  end
  for i = split2, #chromosome2.genes do
    local gene = chromosome2.genes[i]
    if gene ~= nil then table.insert(childChromosome2.genes, copyGene(gene)) end
  end

  --setting maxneurons and mutationRates for both children then return both children
  childChromosome1.maxneuron = math.max(chromosome1.maxneuron, chromosome2.maxneuron)
  childChromosome2.maxneuron = math.max(chromosome1.maxneuron, chromosome2.maxneuron)

  for mutation, rate in pairs(chromosome1.mutationRates) do
    childChromosome1.mutationRates[mutation] = rate
  end
  for mutation, rate in pairs(chromosome1.mutationRates) do
    childChromosome2.mutationRates[mutation] = rate
  end

  return childChromosome1, childChromosome2
end


---------------------------------------------------------------------------------------------------
-- Node mutations applied to chromosomes during mutate functions ----------------------------------
-- Created: March 9, 2016
-- Modified: April 22, 2016
-- Primary author: Christopher Bruns

--Alters a gene/node inside a chromosome
function nodeMutation(chromosome)
  --First check if there are any genes within the chromosome and if not end it
  if #chromosome.genes == 0 then
    return
  end

  --Add one to the maximum number of neurons possible for the chromosome
  chromosome.maxneuron = chromosome.maxneuron + 1

  --Picks a gene randomly from the genes within the chromosome, if it is not enabled the function ends
  local gene = chromosome.genes[math.random(1,#chromosome.genes)]
  if not gene.enabled then
    return
  end

  --Now set it to not enabled since it is being mutated...I THINK...
  gene.enabled = false

  --New gene copy and set the attributes
  local gene1 = copyGene(gene)
  gene1.out = chromosome.maxneuron
  gene1.weight = 1.0
--  gene1.innovation = newInnovation()
  gene1.enabled = true
  --within the table the chromosome has add this new gene
  table.insert(chromosome.genes, gene1)

  --New gene copy and again set some of the attributes
  local gene2 = copyGene(gene)
  gene2.into = chromosome.maxneuron
  gene2.innovation = newInnovation()
  gene2.enabled = true
  --again within the table add the new gene2
  table.insert(chromosome.genes, gene2)

end

--alters a gene/node in the last quarter of the gene table
function tailNodeMutation(chromosome)
  --inside a chromosome, adds new genes to the table w/ .into and .out values
  --First check if there are any genes within the chromosome and if not end it
  if #chromosome.genes == 0 or #chromosome.genes == 1 then
    return
  end

  --Add one to the maximum number of neurons possible for the chromosome
  chromosome.maxneuron = chromosome.maxneuron + 1

  --Picks a gene randomly from the last quarter of genes in the chromosome, if it is not enabled the function ends
  local gene = chromosome.genes[math.random(math.floor(#chromosome.genes / 1.25),#chromosome.genes)]
  if not gene.enabled then
    return
  end

  --Now set it to not enabled since it is being mutated.
  gene.enabled = false
  --New gene copy and set the attributes
  local gene1 = copyGene(gene)
  gene1.out = chromosome.maxneuron
  gene1.weight = 1.0
  gene1.innovation = newInnovation()
  gene1.enabled = true
  --within the table the chromosome has add this new gene
  table.insert(chromosome.genes, gene1)

  --New gene copy and again set some of the attributes
  local gene2 = copyGene(gene)
  gene2.into = chromosome.maxneuron
  gene2.innovation = newInnovation()
  gene2.enabled = true
  --again within the table add the new gene2
  table.insert(chromosome.genes, gene2)
end

--alters two genes/nodes inside a chromosome, adds new genes to the table w/ .into and .out values
function twoNodeMutation(chromosome)
  --First check if there are any genes within the chromosome and if not end it
  if #chromosome.genes == 0 then
    return
  end

  --If there is only one gene call it to do a nodeMutation
  if #chromosome.genes == 1 then
    nodeMutation(chromosome)
    return
  end

  --Add two to the maximum number of neurons possible for the chromosome
  chromosome.maxneuron = chromosome.maxneuron + 2

  --Picks a gene randomly from the genes within the chromosome, if it is not enabled the function ends
  local gene = chromosome.genes[math.random(1,#chromosome.genes)]
  if not gene.enabled then
    return
  end

  --Pick the second gene randomly from the genes within the chromosome
  local secondGene = chromosome.genes[math.random(1,#chromosome.genes)]

  --If it is the same as the first gene do a node mutation
  if secondGene == gene then
    nodeMutation(chromosome)
    return
  end

  --If it is not enabled the function ends
  if not secondGene.enabled then
    return
  end

  --Now set it to not enabled since it is being mutated.
  gene.enabled = false
  --New gene copy and set the attributes
  local gene1 = copyGene(gene)
  gene1.out = chromosome.maxneuron
  gene1.weight = 1.0
  gene1.innovation = newInnovation()
  gene1.enabled = true
  --within the table the chromosome has add this new gene
  table.insert(chromosome.genes, gene1)

  --New gene copy and again set some of the attributes
  local gene2 = copyGene(gene)
  gene2.into = chromosome.maxneuron
  gene2.innovation = newInnovation()
  gene2.enabled = true
  --again within the table add the new gene2
  table.insert(chromosome.genes, gene2)


  --Now set it to not enabled since it is being mutated.
  secondGene.enabled = false
  --New gene copy and set the attributes
  local gene3 = copyGene(secondGene)
  gene3.out = chromosome.maxneuron
  gene3.weight = 1.0
  gene3.innovation = newInnovation()
  gene3.enabled = true
  --within the table the chromosome has add this new gene
  table.insert(chromosome.genes, gene3)

  --New gene copy and again set some of the attributes
  local gene4 = copyGene(secondGene)
  gene4.into = chromosome.maxneuron
  gene4.innovation = newInnovation()
  gene4.enabled = true
  --again within the table add the new gene2
  table.insert(chromosome.genes, gene4)

end

--used to alter the weight of a genes within a chromosome
function stepMutation(chromosome)
  --sets the step to value of step in mutationRates table. Initially set to 0.1
  local step = chromosome.mutationRates["step"]
  --[[Finds the number of genes in the chromosome. For each gene if the random number between 0-1 is less than
  the PerturbChance value (0.9) then the weight of the gene is increased otherwise it is decreased.]]
  for i=1,#chromosome.genes do
    local gene = chromosome.genes[i]
    local random = math.random()
    --[[Will decrease the weight of the gene if it has a lower value (< 1).
    For higher values will increase it some.
    Very little change but really just depends on the previous weight.]]
    if random < PerturbChance then
      gene.weight = gene.weight + random * step*2 - step
    --Will reset the gene weight somewhere between 1.6-2.0
    else
      gene.weight = random*4-2
    end
  end
end


---------------------------------------------------------------------------------------------------
-- Altering links within chromosomes during mutate functions --------------------------------------
-- Created: March 10, 2016
-- Modified: April 20, 2016
-- Primary author: Christopher Bruns


--alters a link between neurons in a chromosome and can be set w/ true or false
function linkMutation(chromosome, forceBias)
  --Get two random neurons from the chromosome given, one of the calls w/ false and the other w/ true
  local neuron1 = randomNeuron(chromosome.genes, false)
  local neuron2 = randomNeuron(chromosome.genes, true)

  --create a newGene
  local newLink = newGene()
  --compare both neurons with the Inputs to see if they are input nodes or not and return if they both are
  if neuron1 <= Inputs and neuron2 <= Inputs then
    return
  end
  --if neuron2 is an input node, then create local variable and set it equal to neuron1. Then set neuron1 equal to
  --the input node neuron2 and set neuron2 to what neuron1 was w/ the variable that was created
  if neuron2 <= Inputs then
     -- Swap output and input
    local temp = neuron1
    neuron1 = neuron2
    neuron2 = temp
  end
  --Now w/ the new gene created (newLink) set the into value to neuron1 and the output value to neuron2
  --Based on the T/F of the forceBias given set it to the Inputs or not
  newLink.into = neuron1
  newLink.out = neuron2
  if forceBias then
    newLink.into = Inputs
  end
  --Check if the genes within the chromosome already has this new gene (newLink) and if it does return
  --Otherwise add the newInnovation and set the weight between 1.6-2.0 and finally add the link to the table
  --of genes within the chromosome
  if containsLink(chromosome.genes, newLink) then
    return
  end
  newLink.innovation = newInnovation()
  newLink.weight = math.random()*4-2

  table.insert(chromosome.genes, newLink)
end

--alters two links between neurons in a chromosome
function doubleLinkMutation(chromosome, forceBias)
  --[[The code for this mutation is very similar to the linkMutation function except it attempts
  to alter two links so comments are less since the implementation is the same except for checking
  if it is able to perform two link mutations and if not it performs only one still.]]
  local neuron1 = randomNeuron(chromosome.genes, false)
  local neuron2 = randomNeuron(chromosome.genes, true)

  local neuron3 = randomNeuron(chromosome.genes, false)
  local neuron4 = randomNeuron(chromosome.genes, true)

  local newLink = newGene()
  local secondNewLink = newGene()

  if (neuron1 <= Inputs and neuron2 <= Inputs) and (neuron3 <= Inputs and neuron4 <= Inputs) then
    return
  end

  if neuron1 <= Inputs and neuron2 <= Inputs then
    if neuron4 <= Inputs then
      local temp = neuron3
      neuron3 = neuron4
      neuron4 = temp
    end

    secondNewLink.into = neuron3
    secondNewLink.out = neuron4
    if forceBias then
      secondNewLink.into = Inputs
    end
    if containsLink(chromosome.genes, secondNewLink) then
      return
    end
    secondNewLink.innovation = newInnovation()
    secondNewLink.weight = math.random()*4-2

    table.insert(chromosome.genes, secondNewLink)
  end

  if neuron3 <= Inputs and neuron4 <= Inputs then
    if neuron2 <= Inputs then
      local temp = neuron1
      neuron1 = neuron2
      neuron2 = temp
    end

    newLink.into = neuron1
    newLink.out = neuron2
    if forceBias then
      newLink.into = Inputs
    end
    if containsLink(chromosome.genes, newLink) then
      return
    end
    newLink.innovation = newInnovation()
    newLink.weight = math.random()*4-2

    table.insert(chromosome.genes, newLink)
  end

  if (neuron1 > Inputs and neuron2 > Inputs) and (neuron3 > Inputs and neuron4 > Inputs) then
     if neuron4 <= Inputs then
      local temp = neuron3
      neuron3 = neuron4
      neuron4 = temp
    end

    secondNewLink.into = neuron3
    secondNewLink.out = neuron4
    if forceBias then
      secondNewLink.into = Inputs
    end
    if containsLink(chromosome.genes, secondNewLink) then
      return
    end
    secondNewLink.innovation = newInnovation()
    secondNewLink.weight = math.random()*4-2

    table.insert(chromosome.genes, secondNewLink)

    if neuron2 <= Inputs then
      local temp = neuron1
      neuron1 = neuron2
      neuron2 = temp
    end

    newLink.into = neuron1
    newLink.out = neuron2
    if forceBias then
      newLink.into = Inputs
    end
    if containsLink(chromosome.genes, newLink) then
      return
    end
    newLink.innovation = newInnovation()
    newLink.weight = math.random()*4-2

    table.insert(chromosome.genes, newLink)
  end
end

--sets genes in a gene eligible or not for mutations
function enableDisableMutations(chromosome, enable)
  --used to put some chromosomes in a table to be chosen from randomly which will change their status to
  --disable if enabled currently
  local candidates = {}
  --Makes a list of chromosomes that are the opposite of the enable variable.
  for _,gene in pairs(chromosome.genes) do
    if gene.enabled == not enable then
      table.insert(candidates, gene)
    end
  end

  if #candidates == 0 then
    return
  end
  --randomly select one of the chromosomes and flip it.
  local gene = candidates[math.random(1,#candidates)]
  gene.enabled = not gene.enabled
end

--checing if a link is within the genes
function containsLink(genes, link)
  --[[for each item in the table of genes provided check if the into and out variables match between the link and the
  gene, if they do return true that it does contain the link.]]
  for i=1,#genes do
    local gene = genes[i]
    if gene.into == link.into and gene.out == link.out then
      return true
    end
  end
end


---------------------------------------------------------------------------------------------------
-- Functions to mutate chromosomes based on set options, random options ---------------------------
-- Created: March 21, 2016
-- Modified: April 22, 2016
-- Primary author: Christopher Bruns

--function called to handle all mutations of a chromosome within
function mutateChromosome(chromosome)
  --[[Takes the chromosome and looks at the mutation table. It randomly (50/50) decides to slightly
  increase or slightly decrease the value of the rate attached with each.]]
  for mutation,rate in pairs(chromosome.mutationRates) do
    if math.random(1,2) == 1 then
      chromosome.mutationRates[mutation] = 0.95*rate
    else
      chromosome.mutationRates[mutation] = 1.05263*rate
    end
  end

  --[[After updating the values in each mutation of the table randomly decided to do a step mutation
  which will update the weight value of each gene in the chromosome]]
  if math.random() < chromosome.mutationRates["connections"] then
    stepMutation(chromosome)
  end

  --[[Likewise will possibly call on the linkMutate which will compare two neurons within the chromosome
  as described in the linkMutate function above.Each time called it decrements the value in the table by 1.]]
  local p = chromosome.mutationRates["link"]
  while p > 0 do
    if math.random() < p then
 
      linkMutation(chromosome, false)
    end
    p = p - 1
  end

  --[[Also will possibly call the linkMutate but with the value of true attached instead of false as with
  the call above with the "link" value. Each time called the value in the bias is decremented by one.]]
  p = chromosome.mutationRates["bias"]
  while p > 0 do
    if math.random() < p then

      linkMutation(chromosome, true)
    end
    p = p - 1
  end

  --[[if the value is greater than 0 it will randomly decide to call nodeMutate which is described above
  each time called it decrements the value within the node mutate table by one and continues]]
  p = chromosome.mutationRates["node"]
  while p > 0 do
    if math.random() < p then
      nodeMutation(chromosome)
    end
    p = p - 1
  end


  p = chromosome.mutationRates["enable"]
  while p > 0 do
    if math.random() < p then
      enableDisableMutations(chromosome, true)
    end
    p = p - 1
  end


  p = chromosome.mutationRates["disable"]
  while p > 0 do
    if math.random() < p then
      enableDisableMutations(chromosome, false)
    end
    p = p - 1
  end
end

--function similar to mutateChromosome but allows for random choices between all mutation options
function randomMutateChromosome(chromosome)
  --[[Takes the chromosome and looks at the mutation table. It randomly (50/50) decides to slightly
  increase or slightly decrease the value of the rate attached with each.]]
  for mutation,rate in pairs(chromosome.mutationRates) do
    if math.random(1,2) == 1 then
      chromosome.mutationRates[mutation] = 0.95*rate
    else
      chromosome.mutationRates[mutation] = 1.05263*rate
    end
  end

  if math.random() < chromosome.mutationRates["connections"] then
    stepMutation(chromosome)
  end

  --[[The next two give the option to randomly run either link or doubleLink mutation]]
  local p = chromosome.mutationRates["link"]
  while p > 0 do
    if math.random() < p then
      if math.random(1,2) == 1 then
        linkMutation(chromosome, false)
      else
        doubleLinkMutation(chromosome, false)
      end
    end
    p = p - 1
  end

  p = chromosome.mutationRates["bias"]
  while p > 0 do
    if math.random() < p then
      if math.random(1,2) == 1 then
        linkMutation(chromosome, true)
      else
        doubleLinkMutation(chromosome, true)
      end
    end
    p = p - 1
  end

  --allows to random choose between the three types of node mutations
  p = chromosome.mutationRates["node"]
  while p > 0 do
    if math.random() < p then
      local randomNum = math.random(1,3)
      if randomNum == 1 then
        nodeMutation(chromosome)
      else if randomNum == 2 then
        twoNodeMutation(chromosome)
      else
        tailNodeMutation(chromosome)
      end
      end
    end
    p = p - 1
  end


  p = chromosome.mutationRates["enable"]
  while p > 0 do
    if math.random() < p then
      enableDisableMutations(chromosome, true)
    end
    p = p - 1
  end


  p = chromosome.mutationRates["disable"]
  while p > 0 do
    if math.random() < p then
      enableDisableMutations(chromosome, false)
    end
    p = p - 1
  end
end

--[[function similar to randomMutate that allows for all types of mutations but is determined based
on the generation number that the population is currently at]]
function generationMutateChromosome(chromosome)
--[[chance to increase or decrease rates but if it increases it increases much more than the other two mutation
functions did]]
  if chromosome.generation < 5 then
    for mutation,rate in pairs(chromosome.mutationRates) do
      if math.random(1,2) == 1 then
        chromosome.mutationRates[mutation] = 0.95*rate
      else
        chromosome.mutationRates[mutation] = 1.25263*rate
      end
    end
  else
    for mutation,rate in pairs(chromosome.mutationRates) do
      if math.random(1,2) == 1 then
        chromosome.mutationRates[mutation] = 0.95*rate
      else
        chromosome.mutationRates[mutation] = 1.05263*rate
      end
    end
  end



  if math.random() < chromosome.mutationRates["connections"] then
    stepMutation(chromosome)
  end


  local p = chromosome.mutationRates["link"]
  while p > 0 do
    if math.random() < p then
      linkMutation(chromosome, false)
    end
    p = p - 1
  end


  p = chromosome.mutationRates["bias"]
  while p > 0 do
    if math.random() < p then
      linkMutation(chromosome, true)
    end
    p = p - 1
  end

  --node mutations change based on the generation number
  p = chromosome.mutationRates["node"]
  while p > 0 do
    if math.random() < p then
      if chromosome.generation > 5 then tailNodeMutation(chromosome)
      else
        local ran = math.random(1,2)
        if ran == 1 then
          nodeMutation(chromosome)
        else twoNodeMutation(chromosome)
        end
      end
    end
    p = p - 1
  end


  p = chromosome.mutationRates["enable"]
  while p > 0 do
    if math.random() < p then
      enableDisableMutations(chromosome, true)
    end
    p = p - 1
  end


  p = chromosome.mutationRates["disable"]
  while p > 0 do
    if math.random() < p then
      enableDisableMutations(chromosome, false)
    end
    p = p - 1
  end
end


---------------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------
-- Interaction with BizHawk -----------------------------------------------------------------------
-- Last modified: March 22, 2016
-- Primary author: Jacob Espenscheid

function getSprites()
	local sprites = {}
	for slot=0,11 do
		local status = memory.readbyte(0x14C8+slot)
		if status ~= 0 then
			spritex = memory.readbyte(0xE4+slot) + memory.readbyte(0x14E0+slot)*256
			spritey = memory.readbyte(0xD8+slot) + memory.readbyte(0x14D4+slot)*256
			sprites[#sprites+1] = {["x"]=spritex, ["y"]=spritey}
		end
	end

	return sprites
end


function getExtendedSprites()
	local extended = {}
	for slot=0,11 do
		local number = memory.readbyte(0x170B+slot)
		if number ~= 0 then
			spritex = memory.readbyte(0x171F+slot) + memory.readbyte(0x1733+slot)*256
			spritey = memory.readbyte(0x1715+slot) + memory.readbyte(0x1729+slot)*256
			extended[#extended+1] = {["x"]=spritex, ["y"]=spritey}
		end
	end

	return extended
end

--location of mario in the environment
function getTile(dx, dy)
	x = math.floor((marioX+dx+8)/16)
	y = math.floor((marioY+dy)/16)

	return memory.readbyte(0x1C800 + math.floor(x/0x10)*0x1B0 + y*0x10 + x%0x10)
end

--locations and inputs based on his location as well
function getInputs()
	getPositions()

	sprites = getSprites()
	extended = getExtendedSprites()

	local inputs = {}

	for dy=-BoxRadius*16,BoxRadius*16,16 do
		for dx=-BoxRadius*16,BoxRadius*16,16 do
			inputs[#inputs+1] = 0

      --get the byte value for each "tile" in the frame
			tile = getTile(dx, dy)
			if tile == 1 and marioY+dy < 0x1B0 then
				inputs[#inputs] = 1
			end

      --This sets the input value to -1 if Mario is near a sprite
			for i = 1,#sprites do
				distx = math.abs(sprites[i]["x"] - (marioX+dx))
				disty = math.abs(sprites[i]["y"] - (marioY+dy))
				if distx <= 8 and disty <= 8 then
					inputs[#inputs] = -1
				end
			end

      --same as above
			for i = 1,#extended do
				distx = math.abs(extended[i]["x"] - (marioX+dx))
				disty = math.abs(extended[i]["y"] - (marioY+dy))
				if distx < 8 and disty < 8 then
					inputs[#inputs] = -1
				end
			end
		end
	end

  return inputs
end

function clearJoypad()
  controller = {}
  for b = 1,#ButtonNames do
    controller["P1 " .. ButtonNames[b]] = false
  end
  joypad.set(controller)
end

--getting current position of Mario in the emulator
function getPositions()
	marioX = memory.read_s16_le(0x94)
	marioY = memory.read_s16_le(0x96)

	local layer1x = memory.read_s16_le(0x1A);
	local layer1y = memory.read_s16_le(0x1C);

	screenX = marioX-layer1x
	screenY = marioY-layer1y
end

---------------------------------------------------------------------------------------------------
-- Crossover in Breed Child -----------------------------------------------------------------------
-- Last modified: April 13, 2016
-- Primary author: Christopher Bruns

--called to breed new children with the crossover technique then they are mutated before being returned
function breedChild(species)
	local child = {}
	--get to chromosomes from the species and crossover with them
	if math.random() < CrossoverChance then
		g1 = species.chromosomes[math.random(1, #species.chromosomes)]
		g2 = species.chromosomes[math.random(1, #species.chromosomes)]
		child = crossoverMutation(g1, g2)
	else
		g = species.chromosomes[math.random(1, #species.chromosomes)]
		child = copyChromosome(g)
	end

	mutateChromosome(child)

	return child
end


---------------------------------------------------------------------------------------------------
-- Alternate Crossover Options in Breed Child -----------------------------------------------------
-- Last modified: April 25, 2016
-- Primary author: Christopher Bruns

--alternate breen child that picks between the three crossover techniques randomly and uses generationMutate function
function altBreedChild(species)
  local child = {}
  local child2 = {}

  if math.random() < CrossoverChance then
    local ranMath = math.random(1,3)
    if ranMath == 1 then
      g1 = species.chromosomes[math.random(1, #species.chromosomes)]
      g2 = species.chromosomes[math.random(1, #species.chromosomes)]
      child = crossoverMutation(g1, g2)
    else if ranMath == 2 then
      g1 = species.chromosomes[math.random(1, #species.chromosomes)]
      g2 = species.chromosomes[math.random(1, #species.chromosomes)]
      child = onePointCrossover(g1, g2)
    else
      g1 = species.chromosomes[math.random(1, #species.chromosomes)]
      g2 = species.chromosomes[math.random(1, #species.chromosomes)]
      child = twoPointCrossover(g1, g2)
    end
  end
  else
    c = species.chromosomes[math.random(1, #species.chromosomes)]
    child = copyChromosome(c)
  end

  generationMutateChromosome(child)

  return child
end

---------------------------------------------------------------------------------------------------
-- Neural network evaluation ----------------------------------------------------------------------
-- Last modified: March 29, 2016
-- Primary author: Jack Van Gent

--function to evaluate the network
function evaluateNetwork(network, inputs)
	table.insert(inputs, 1)
	if #inputs ~= Inputs then
		console.writeline("Incorrect number of neural network inputs.")
		return {}
	end
	--Input = 170
	for i=1,Inputs do
		network.neurons[i].value = inputs[i]
	end
	--loop through neurons
	for _,neuron in pairs(network.neurons) do
		local sum = 0
		--loop through the incoming neurons
		for j = 1,#neuron.incomingValues do
			--get incoming gene
			local incoming = neuron.incomingValues[j]
			-- get other neuron
			local other = network.neurons[incoming.into]
			--sum is a running total of weight * value
			sum = sum + incoming.weight * other.value
		end

		if #neuron.incomingValues > 0 then
			--set neuron's value to sigmoid of the sum
			neuron.value = sigmoid(sum)
		end
	end
	--creates a new list of outputs to the controller
	local outputs = {}
	for o=1,Outputs do
		--creates new button based on buttonNames
		local button = "P1 " .. ButtonNames[o]
		if network.neurons[MaxNeurons+o].value > 0 then
			outputs[button] = true
		else
			outputs[button] = false
		end
	end

	return outputs
end

--used to evalute the current network
function evaluateCurrent()
  --this function traverses a network and produces a move from it for Mario to make
  local species = population.species[population.currentSpecies]
  local chromosome = species.chromosomes[population.currentChromosome]

  inputs = getInputs()
  --gets a new set of outputs for controller
  controller = evaluateNetwork(chromosome.network, inputs)
  --if both left and right are pressed unpress them
  if controller["P1 Left"] and controller["P1 Right"] then
    controller["P1 Left"] = false
    controller["P1 Right"] = false
  end
  --if both up and down are pressed unpress them
  if controller["P1 Up"] and controller["P1 Down"] then
    controller["P1 Up"] = false
    controller["P1 Down"] = false
  end

  joypad.set(controller)
end

function sigmoid(x)
	return 2/(1+math.exp(-4.9*x))-1
end


---------------------------------------------------------------------------------------------------
-- Preparation for a new generation ---------------------------------------------------------------
-- Last modified: March 22, 2016
-- Primary author: Jacob Espenscheid

--removing species that have not been changing or mutating much and their fitness is also poor
function removeStaleSpecies()
	--list for survivors
	local survived = {}
	--loop through population species
	for s = 1,#population.species do
		local species = population.species[s]
		--sort species chromosomes
		table.sort(species.chromosomes, function (a,b)
			return (a.fitness > b.fitness)
		end)
		--if the first chromosome in the list has a higher fitness
		if species.chromosomes[1].fitness > species.topFitness then
			species.topFitness = species.chromosomes[1].fitness
			--resets staleness to zero because it did better
			species.staleness = 0
		else
			--add one because it didn't do better
			species.staleness = species.staleness + 1
		end
		--if staleness < 15 or species top fitness is >= the best keep it
		--if the species hasn't changed in 15 rounds and isnt the best it's culled
		if species.staleness < StaleSpecies or species.topFitness >= population.maxFitness then
			table.insert(survived, species)
		end
	end
	--make population.species = survived
	population.species = survived
end

--removing weak species based on the average fitness
function removeWeakSpecies()
	--list for survived
	local survived = {}
	--totalAverageFitness returns the sum of each species avgerage fitness
	local sum = totalAverageFitness()
	--for each species
	for s = 1,#population.species do
		local species = population.species[s]
		--we set the population as a global in the beginning
		--if they pass this if statement they live
		breed = math.floor(species.averageFitness / sum * MaxPopulation)
		if breed >= 1 then
			table.insert(survived, species)
		end
	end
	--make the list of survived the population.species
	population.species = survived
end

--sets up a new generation
function newGeneration()
	-- Cull the bottom half of each species
	cullSpecies(false)
	--ranks all chromosomes by fitness
	rankGlobally()
	--culls species that haven't changed in a long time and aren't the best
	removeStaleSpecies()
	--rank again after species could have been removed
	rankGlobally()
	--loop through species
	for s = 1,#population.species do
		local species = population.species[s]
		--calculate average fitness for each species
		calculateAverageFitness(species)
	end
	--uses average fitness to get rid of weak species
	removeWeakSpecies()
	--sum is the all the species average fitness totaled
	local sum = totalAverageFitness()
	--new list children for the next species
	local children = {}
	for s = 1,#population.species do
		local species = population.species[s]
		breed = math.floor(species.averageFitness / sum * MaxPopulation) - 1
		for i=1,breed do
			table.insert(children, breedChild(species))
		end
	end
	cullSpecies(true) -- Cull all but the top member of each species
	while #children + #population.species < MaxPopulation do
		local species = population.species[math.random(1, #population.species)]
		table.insert(children, breedChild(species))
	end
	for c=1,#children do
		local child = children[c]
		addToSpecies(child)
	end

	population.generation = population.generation + 1
  for i = 1, #population.species do
    local species = population.species[i]
    for x = 1, #species.chromosomes do
      local chromosome = species.chromosomes[x]
      chromosome.generation = population.generation
    end
  end
	writeFile("backup." .. population.generation .. "." .. DatafileName)
end


function fitnessAlreadyMeasured()
	local species = population.species[population.currentSpecies]
	local chromosome = species.chromosomes[population.currentChromosome]

	return chromosome.fitness ~= 0
end

--grabs the next chromosome to be used
function nextChromosome()
  --grab next chromosome
	population.currentChromosome = population.currentChromosome + 1

  --if we're at the last chromosome in the species, start a new one and grab first chromosome
	if population.currentChromosome > #population.species[population.currentSpecies].chromosomes then
		population.currentChromosome = 1
		population.currentSpecies = population.currentSpecies+1
		if population.currentSpecies > #population.species then
			newGeneration()
			population.currentSpecies = 1
		end
	end
end

--ranking all within the population in order of fitness
function rankGlobally()
	--list for global
	local global = {}
	--loop through the species
	for s = 1,#population.species do
		local species = population.species[s]
		--add each chromosome into the global list
		for g = 1,#species.chromosomes do
			table.insert(global, species.chromosomes[g])
		end
	end
	--sort global based on fitness
	table.sort(global, function (a,b)
		return (a.fitness < b.fitness)
	end)
	--for each chromosome in global
	for g=1,#global do
		--set that chromosomes rank to its position in the list
		global[g].globalRank = g
	end
end

--finds the average fitness of a species
function calculateAverageFitness(species)
	local total = 0

	for g=1,#species.chromosomes do
		local chromosome = species.chromosomes[g]
		total = total + chromosome.globalRank
	end

	species.averageFitness = total / #species.chromosomes
end

--finds the average fitness of the entire population
function totalAverageFitness()
	local total = 0
	for s = 1,#population.species do
		local species = population.species[s]
		total = total + species.averageFitness
	end

	return total
end

--if true is sent in removes all but the highest rated chromosome in a species
function cullSpecies(cutToOne)
	--Loop through population species
	for s = 1,#population.species do
		--select species
		local species = population.species[s]

		--sort chromosomes by fitness
		table.sort(species.chromosomes, function (a,b)
			return (a.fitness > b.fitness)
		end)
		--remaining is how many are left after cull
		local remaining = math.ceil(#species.chromosomes/2)
		--cutToOne is either true or false
		if cutToOne then
			remaining = 1
		end
		--loop through and remove the last chromosome in the table until you hit remainder
		while #species.chromosomes > remaining do
			table.remove(species.chromosomes)
		end
	end
end


---------------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------
-- File IO for data backups -----------------------------------------------------------------------
-- Last modified: April 16, 2016
-- Primary author: Jacob Espenscheid

--[[writing a file containing information on a given population to be reviewed or loaded back later
if the user wishes to continue from this point in a future test or wishes to view information on the
population that was saved]]
function writeFile(backupFileName)
  local file = io.open(backupFileName, "w")
	file:write(population.generation .. "\n")
	file:write(population.maxFitness .. "\n")
	file:write(#population.species .. "\n")
  for n,species in pairs(population.species) do
		file:write(species.topFitness .. "\n")
		file:write(species.staleness .. "\n")
		file:write(#species.chromosomes .. "\n")
		for m,chromosome in pairs(species.chromosomes) do
			file:write(chromosome.fitness .. "\n")
			file:write(chromosome.maxneuron .. "\n")
			file:write(chromosome.generation .. "\n")
			for mutation,rate in pairs(chromosome.mutationRates) do
				file:write(mutation .. "\n")
				file:write(rate .. "\n")
			end
			file:write("done\n")

			file:write(#chromosome.genes .. "\n")
			for l,gene in pairs(chromosome.genes) do
				file:write(gene.into .. " ")
				file:write(gene.out .. " ")
				file:write(gene.weight .. " ")
				file:write(gene.innovation .. " ")
				if(gene.enabled) then
					file:write("1\n")
				else
					file:write("0\n")
				end
			end
		end
  end
  file:close()
end

--[[function to read and and load a file containing a previously saved population to continue where
a previous run left off]]
function loadFile(filename)
	if filename == "" then
		return
	end
  local file = io.open(filename, "r")
	population = newPopulation()
	population.generation = file:read("*number")
	population.maxFitness = file:read("*number")
	--forms.settext(maxFitnessLabel, "Max Fitness: " .. math.floor(pool.maxFitness))
	local numSpecies = file:read("*number")
  for s=1,numSpecies do
		local species = newSpecies()
		table.insert(population.species, species)
		species.topFitness = file:read("*number")
		species.staleness = file:read("*number")
		local numChromosomes = file:read("*number")
		for g=1,numChromosomes do
			local chromosome = newChromosome()
			table.insert(species.chromosomes, chromosome)
			chromosome.fitness = file:read("*number")
			chromosome.maxneuron = file:read("*number")
			chromosome.generation = file:read("*number")
			local line = file:read("*line")
			while line ~= "done" do
				chromosome.mutationRates[line] = file:read("*number")
				line = file:read("*line")
			end
			local numGenes = file:read("*number")
			for n=1,numGenes do
				local gene = newGene()
				table.insert(chromosome.genes, gene)
				local enabled
				gene.into, gene.out, gene.weight, gene.innovation, enabled = file:read("*number", "*number", "*number", "*number", "*number")
				if enabled == 0 then
					gene.enabled = false
				else
					gene.enabled = true
				end

			end
		end
	end
  file:close()

	while fitnessAlreadyMeasured() do
		nextChromosome()
	end
	initializeRun()
	population.currentFrame = population.currentFrame + 1
end

---------------------------------------------------------------------------------------------------
-- File IO for test data --------------------------------------------------------------------------
-- Last modified: April 2, 2016
-- Primary author: Jack Van Gent

--formating the time for when writing to the file
function formatTime(time)
  local hours = math.floor(time / 3600)
  local leftoverSeconds = (time%3600)
  local minutes = math.floor(leftoverSeconds / 60)
  local seconds = math.floor(leftoverSeconds%60)
  return hours, minutes, seconds
end

--writing information to file based on milestone performance
function writeMilestone(milestone, fitness)
  currentTime = os.clock() - StartTime
  hours, minutes, seconds = formatTime(currentTime)
  openFile()
  io.write(milestone, ": ", hours, "h, ", minutes, "m, ", seconds, "s. FITNESS: ", fitness, "\n")
  File:close()
end

--[[formating the output that is written into files and checking the information with
the milestones it reached and the time it took to reach those milestones, specifically
when it has completed a level successfully]]
function checkMilestone(fitness, won)
  if testDataOn then
    if won then
      currentTime = os.clock() - StartTime
      hours, minutes, seconds = formatTime(currentTime)
      openFile()
      io.write("COMPLETED: ", hours, "h, ", minutes, "m, ", seconds, "s. FITNESS: ", fitness, "\n")
      io.write("END\n")
      File:close()
    else
      currentTime = os.clock() - StartTime
      for i=1,#Milestones do
        if Milestones[i].reached ~= true then
          if fitness >= Milestones[i].number then
            Milestones[i].reached = true
            writeMilestone(Milestones[i].number, fitness)
          end
        end
      end
    end
  else
    return
  end
end

function openFile()
  File = io.open(DatafileName, "a+")
  io.output(File)
end

--create milestones to keep track of progress as Mario reaches different points in the level based on fitness score
function generateMilestones()
  Milestones = {}
  for i=100,5000,100 do
    entry = {number=i, reached=false}
    table.insert(Milestones, entry)
  end
end

--set the IO for files created when testing different types of mutation combinations
function initializeIO()
  generateMilestones()
  openFile()
  io.write(Description, "\n")
  io.write("START\n")
  StartTime = os.clock() --returns the number of seconds of CPU time for the program (in seconds)
  File:close()
end


---------------------------------------------------------------------------------------------------
---------------------------------------------------------------------------------------------------
-- Data initialization ----------------------------------------------------------------------------
-- Last modified: March 7, 2016
-- Primary author: Jacob Espenscheid

--initializing the running of the system and creating a network
function initializeRun()
  savestate.load(Filename);
  rightmost = 0
  population.currentFrame = 0
  timeout = TimeoutConstant
  clearJoypad()

  local species = population.species[population.currentSpecies]
  local chromosome = species.chromosomes[population.currentChromosome]
  createNetwork(chromosome)
  evaluateCurrent()
end

--set up an initial population
function initializePopulation()
  population = newPopulation()

  --create initial population of chromosomes, empty for now
  for i=1,MaxPopulation do
    tempChromosome = basicChromosome()
    addToSpecies(tempChromosome)
  end

  initializeRun()
end


---------------------------------------------------------------------------------------------------
-- Misc functions ---------------------------------------------------------------------------------
-- Last modified: April 22, 2016
-- Primary author: All

--used to add a chromosome child to a species
function addToSpecies(child)
  local foundSpecies = false
  for s=1,#population.species do
    local species = population.species[s]
    if not foundSpecies and sameSpecies(child, species.chromosomes[1]) then
      table.insert(species.chromosomes, child)
      foundSpecies = true
    end
  end

  if not foundSpecies then
    local childSpecies = newSpecies()
    table.insert(childSpecies.chromosomes, child)
    table.insert(population.species, childSpecies)
  end
end

--taking chromosomes from same species and changing weights and comparing weights in genes from the chromosomes
function sameSpecies(chromosome1, chromosome2)
  local dd = DeltaDisjoint*disjoint(chromosome1.genes, chromosome2.genes)
  local dw = DeltaWeights*weights(chromosome1.genes, chromosome2.genes)
  return dd + dw < DeltaThreshold
end

function disjoint(genes1, genes2)
  local i1 = {}
  for i = 1,#genes1 do
    local gene = genes1[i]
    i1[gene.innovation] = true
  end

  local i2 = {}
  for i = 1,#genes2 do
    local gene = genes2[i]
    i2[gene.innovation] = true
  end

  local disjointGenes = 0
  for i = 1,#genes1 do
    local gene = genes1[i]
    if not i2[gene.innovation] then
      disjointGenes = disjointGenes+1
    end
  end

  for i = 1,#genes2 do
    local gene = genes2[i]
    if not i1[gene.innovation] then
      disjointGenes = disjointGenes+1
    end
  end

  local n = math.max(#genes1, #genes2)

  return disjointGenes / n
end

--comparing and setting weights of genes
function weights(genes1, genes2)
  local i2 = {}
  for i = 1,#genes2 do
    local gene = genes2[i]
    i2[gene.innovation] = gene
  end

  local sum = 0
  local coincident = 0
  for i = 1,#genes1 do
    local gene = genes1[i]
    if i2[gene.innovation] ~= nil then
      local gene2 = i2[gene.innovation]
      sum = sum + math.abs(gene.weight - gene2.weight)
      coincident = coincident + 1
    end
  end

  return sum / coincident
end


---------------------------------------------------------------------------------------------------
-- Startup activities
if testDataOn then
  initializeIO()
end

if population == nil then
  initializePopulation()
end

loadFile(loadFileName)


-- MAIN LOOP --------------------------------------------------------------------------------------
-- Last modified: April 24, 2016
-- Primary author: Jack Van Gent

while true do

  local species = population.species[population.currentSpecies]
  local chromosome = species.chromosomes[population.currentChromosome]

  --evaluate current every 5 frames
  if population.currentFrame%5 == 0 then
		evaluateCurrent()
	end

  --resends the move produced by evaluateCurrent to bizhawk
  joypad.set(controller)

  --we need to know mario's position for the timeout
  --keep resetting timeout if Mario keeps moving right
  getPositions()
	if marioX > rightmost then
		rightmost = marioX
		timeout = TimeoutConstant
	end

  timeout = timeout - 1

  --code below happens when level is beaten or timeout happens (or if he dies)
  local timeoutBonus = population.currentFrame / 4
	if timeout + timeoutBonus <= 0 then

    --this is the fitness function, distance/time
		local fitness = rightmost - population.currentFrame / 2

    --did he beat the level?
		if rightmost > 4816 then
			fitness = fitness + 1000
      checkMilestone(fitness, true)
    else
      checkMilestone(fitness, false)
		end
		if fitness == 0 then
			fitness = -1
		end
		chromosome.fitness = fitness

    --happens on new max fitness
		if fitness > population.maxFitness then
			population.maxFitness = fitness
		end

		console.writeline("Gen " .. population.generation .. " species " .. population.currentSpecies .. " chromosome " .. population.currentChromosome .. " fitness: " .. fitness)

		while fitnessAlreadyMeasured() do
			nextChromosome()
		end
		initializeRun()
	end

  population.currentFrame = population.currentFrame + 1
  emu.frameadvance();
end


---------------------------------------------------------------------------------------------------

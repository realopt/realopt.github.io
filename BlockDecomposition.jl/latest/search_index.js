var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#Home-1",
    "page": "Home",
    "title": "Home",
    "category": "section",
    "text": "BlockDecomposition.jl package provides features to take advantage of the shape of block structured problems.The user must be familiar with the syntax of JuMP, which is described in its documentation.note: Note\nThis package is under development. Although it is a JuMP extension, it is not written nor maintained by the primary developers of JuMP. Therefore, do not expect high reactiveness on the issues."
},

{
    "location": "index.html#Manual-Outline-1",
    "page": "Home",
    "title": "Manual Outline",
    "category": "section",
    "text": "Pages = [\n    \"index.md\",\n    \"installation.md\",\n    \"introduction.md\",\n    \"basic.md\",\n    \"callbacks.md\",\n    \"advanced.md\",\n    \"BlockSolverInterface.md\"\n]\nDepth = 1"
},

{
    "location": "installation.html#",
    "page": "Installation",
    "title": "Installation",
    "category": "page",
    "text": ""
},

{
    "location": "installation.html#Installation-Guide-1",
    "page": "Installation",
    "title": "Installation Guide",
    "category": "section",
    "text": "BlockDecomposition requires Julia, JuMP.jl and MathProgBase.jl"
},

{
    "location": "installation.html#Getting-BlockDecomposition.jl-1",
    "page": "Installation",
    "title": "Getting BlockDecomposition.jl",
    "category": "section",
    "text": "BlockDecomposition.jl can be installed using the package manager of Julia. To install it, runjulia> Pkg.clone(\"git@github.com:realopt/BlockDecomposition.jl.git\")This command will, recursively, install BlockDecomposition.jl and its dependencies.To start using BlockDecomposition.jl, it should be imported together with JuMP into the local scopeusing JuMP, BlockDecomposition"
},

{
    "location": "introduction.html#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "introduction.html#Introduction-1",
    "page": "Introduction",
    "title": "Introduction",
    "category": "section",
    "text": "BlockDecomposition.jl is a package providing features to take advantage of the shape of block structured problem; in other words, problems on which Dantzig-Wolfe decomposition or Benders decomposition can be applied.note: Note\nThis package provides:     1. a JuMP modeling layer for describing Benders decomposition     2. a JuMP modeling layer for describing Dantzig-Wolfe decomposition     3. an interface to write customized oracles for (Benders/Dantzig-Wolfe) subproblemsOnly the first feature is supported by a solver for now (CPLEX 12.7)"
},

{
    "location": "introduction.html#Benders-decomposition-1",
    "page": "Introduction",
    "title": "Benders decomposition",
    "category": "section",
    "text": "Consider the following problembeginarrayc c c c c c c c\ntextmin        sum_alpha=1^h f_alpha y_alpha       +  sum_alpha=1^i f_1alpha x_1alpha  +  sum_alpha=1^j f_2alpha x_2alpha  +  sum_alpha=1^k f_3alpha x_3alpha        \ntextc_11   sum_alpha=1^h a_11alpha y_alpha  +  sum_alpha=1^i b_1alpha x_1alpha                                                                                                  geq beta_11 \nvdots                                                          vdots                                                                                                                                             \ntextc_1n  sum_alpha=1^h a_1nalpha y_alpha  +  sum_alpha=1^i b_nalpha x_1alpha                                                                                                  geq beta_1n \ntextc_21  sum_alpha=1^h a_21alpha y_alpha  +                                                  sum_alpha=1^j c_1alpha x_2alpha                                                  geq beta_21 \nvdots                                                                                                          vdots                                                                                             \ntextc_2p  sum_alpha=1^h a_2palpha y_alpha  +                                                  sum_alpha=1^j c_palpha x_2alpha                                                  geq beta_2p \ntextc_31  sum_alpha=1^h a_31alpha y_alpha  +                                                                                                  sum_alpha=1^k d_1alpha x_3alpha  geq beta_31 \nvdots                                                                                                                                                          vdots                                             \ntextc_3q  sum_alpha=1^h a_3qalpha y_alpha  +                                                                                                  sum_alpha=1^k d_qalpha x_3alpha  geq beta_3q \nendarrayThe coefficients matrix has the following block structure :(Image: BendersImg)Hence, we can apply Benders decomposition on this problem.We partition the variables.Variables y_alpha, alphain 1 ldots h are in the master.\nVariables x_1alpha, alpha in 1 ldots i are in the first subproblem.\nVariables x_2alpha, alpha in 1 ldots j are in the second subproblem.\nVariables x_3alpha, alpha in 1 ldots k are in the third subproblem.Assume that m is the compact formulation written with JuMP.The user must create a function to describe this decomposition. Such function could be:function b_decomp(var_name, var_id)\n    if var_name == :y\n        return (:B_MASTER, 0)\n    else\n        return (:B_SP, var_id[1])\n    end\nend\nadd_Benders_decomposition(m, b_decomp)"
},

{
    "location": "introduction.html#Dantzig-Wolfe-decomposition-1",
    "page": "Introduction",
    "title": "Dantzig-Wolfe decomposition",
    "category": "section",
    "text": "Consider the following problem :beginarrayc c c c c c c\ntextminimize   sum_alpha=1^i f_alpha x_alpha      +  sum_alpha=1^j f_i+alpha x_i+alpha      +  sum_alpha=1^k f_i+j+alpha x_i+j+alpha                      \ntextmc_1    sum_alpha=1^i a_1alpha x_alpha  +  sum_alpha=1^j a^1_i+alpha x_i+alpha  +  sum_alpha=1^k a_1i+j+alpha x_i+j+alpha  geq  beta_01 \nvdots            vdots                                      vdots                                              vdots                                                               \ntextmc_m    sum_alpha=1^i a_malpha x_alpha  +  sum_alpha=1^j a^m_i+alpha x_i+alpha  +  sum_alpha=1^k a_mi+j+alpha x_i+j+alpha  geq  beta_0m \ntextsc_11  sum_alpha=1^i b_1alpha x_alpha                                                                                                            geq  beta_11 \nvdots            vdots                                                                                                                                                               \ntextsc_1n  sum_alpha=1^i b_nalpha x_alpha                                                                                                            geq  beta_1n \ntextsc_21                                              sum_alpha=1^j c_1i+alpha x_i+alpha                                                          geq  beta_21 \nvdots                                                        vdots                                                                                                                   \ntextsc_2p                                              sum_alpha=1^j c_pi+alpha x_i+alpha                                                          geq  beta_2p \ntextsc_31                                                                                                  sum_alpha=1^k d_1i+j+alpha x_i+j+alpha  geq  beta_31 \nvdots                                                                                                            vdots                                                               \ntextsc_3q                                                                                                  sum_alpha=1^k d_qi+j+alpha x_i+j+alpha  geq  beta_3q \nendarrayThe coefficients matrix has the following block structure :(Image: DWImg)Hence, we can apply Dantzig-Wolfe decomposition on this problem.We partition the constraints.Constraints textmc_1 to textmc_m are in the master.\nConstraints textsc_11 to textsc_1n are in the first subproblem.\nConstraints textsc_21 to textsc_2p are in the second subproblem.\nConstraints textsc_31 to textsc_3q are in the third subproblem.As for Benders decomposition, the user must to create function to describe this decomposition. Such function could be: ::function dw_decomp(constr_name, constr_id)\n    if constr_name == :mc\n        return (:DW_MASTER, 0)\n    else\n        return (:DW_SP, constr_id[1])\n    end\nend\nadd_Dantzig_Wolfe_decomposition(m, dw_decomp)"
},

{
    "location": "basic.html#",
    "page": "Basic",
    "title": "Basic",
    "category": "page",
    "text": ""
},

{
    "location": "basic.html#Basic-usage-1",
    "page": "Basic",
    "title": "Basic usage",
    "category": "section",
    "text": "This quick start guide introduces features of BlockDecomposition.jl package."
},

{
    "location": "basic.html#BlockModel-instantiation-1",
    "page": "Basic",
    "title": "BlockModel instantiation",
    "category": "section",
    "text": "A BlockDecomposition model can be instantiated as  gap = BlockModel(solver = decomp_solver)The instantiation is similar to the one of JuMP.Model. However decomp_solver must be a MIP solver of JuMP models that additionally supports Benders and/or Dantzig-Wolfe decomposition."
},

{
    "location": "basic.html#Write-the-model-1",
    "page": "Basic",
    "title": "Write the model",
    "category": "section",
    "text": "The model is written as a JuMP model. If you are not familiar with JuMP syntax, you may want to check its documentation.The following example is the capacitated facility location problem. Consider a set of potential facility sites Facilities = 1:F where a facility can be opened and a set of customers Customers = 1:C that must be serviced. Assign a customer c to a facility f has a cost DistanceCosts[c, f]. Moreover, opening a facility has a cost Fixedcosts[f] and each facility has a capacity Capacities[f]. All customers must be assigned to a facility.fl = BlockModel(solver = decomp_solver)\n\n@variable(fl, 0 <= x[c in Customers, f in Factories] <= 1 )\n@variable(fl, y[f in Facilities], Bin)\n\n@constraint(fl, cov[c in Customers],\n            sum( x[c, f] for j in Facilities ) >= 1)\n\n@constraint(fl, knp[f in Facilities],\n            sum( x[c, f] for c in Customers ) <= y[f] * Capacities[f])\n\n@objective(fl, Min,\n            sum( DistanceCosts[c, f] * x[c, f] for f in Facilities, c in Customers)\n            + sum( Fixedcosts[f] * y[f] for f in Facilities) )"
},

{
    "location": "basic.html#BlockDecomposition.add_Benders_decomposition",
    "page": "Basic",
    "title": "BlockDecomposition.add_Benders_decomposition",
    "category": "function",
    "text": "add_Benders_decomposition(model::JuMP.Model, B_decomp::Function)\n\nassigns the decomposition function B_decomp to the model.\n\n\n\n"
},

{
    "location": "basic.html#Decomposition-1",
    "page": "Basic",
    "title": "Decomposition",
    "category": "section",
    "text": "The decomposition is described with a function that takes two arguments. This function is call by BlockDecomposition to build the decomposition data. For Benders decomposition, the arguments will be varname the name of the variable and varid the index of the variable.function B_decomp(varname::Symbol, varid::Tuple) :: Tuple{Symbol, Union{Int, Tuple}}The function returns a Tuple that contains a Symbol and a Union{Int, Tuple}. The Symbol is the type of problem to which the variable belongs. It may be :B_MASTER or :B_SP. The Union{Int, Tuple} is the index of this problem.This function should be attached to the BlockModel usingadd_Benders_decompositionNow, we can write the decomposition function of our two ewamples. For the Capacitated Facility Location problem, we want to put variables y in the master and variables x in the unique subproblem. It can be writen asfunction benders_fct(varname::Symbol, varid::Tuple) :: Tuple{Symbol, Union{Int, Tuple}}\n    if varname == :x              # variables x will be assigned to the\n        return (:B_SP, 0)        # subproblem that has the index 0\n    else                          # variables y will be assigned to the\n        return (:B_MASTER, 0)    # master that has the index 0\n    end\nend\nadd_Benders_decomposition(fl, benders_fct)Notice that even if there is only one master problem, he must have an index."
},

{
    "location": "callbacks.html#",
    "page": "Callbacks",
    "title": "Callbacks",
    "category": "page",
    "text": ""
},

{
    "location": "callbacks.html#Callbacks-1",
    "page": "Callbacks",
    "title": "Callbacks",
    "category": "section",
    "text": "note: Note\nNo solver over Julia supports this feature yet.Oracles are customized solvers that can be used to solve efficiently subproblems. We introduce them using the example of Generalized Assignment Problem."
},

{
    "location": "callbacks.html#Introduction-1",
    "page": "Callbacks",
    "title": "Introduction",
    "category": "section",
    "text": "Consider a set of machines Machines = 1:M and a set of jobs Jobs = 1:J. A machine m has a resource capacity Capacity[m]. When we assign a job j to a machine m, the job has a cost Cost[m,j] and consumes Weight[m,j] resources of the machine m. The goal is to minimize the jobs cost sum by assigning each job to a machine while not exceeding the capacity of each machine. The model isgap = BlockModel(solver = solver)\n\n@variable(gap, x[m in Machines, j in Jobs], Bin)\n\n@constraint(gap, cov[0, j in Jobs],\n               sum(x[m,j] for m in Machines) >= 1)\n\n@constraint(gap, knp[m in Machines],\n               sum(Weight[m,j]*x[m,j] for j in Jobs) <= Capacity[m])\n\n@objective(gap, Min,\n               sum(Cost[m,j]*x[m,j] for m in Machines, j in Jobs))\n\nfunction dw_fct(cstrname::Symbol, cstrid::Tuple) :: Tuple{Symbol, Tuple}\n    if cstrname == :cov           # cov constraints will be assigned in the\n        return (:DW_MASTER, (0,)) # master that has the index 0\n    else                          # others constraints will be assigned in a\n        return (:DW_SP, cstrid)   # subproblem with same index as the constraint\n    end\nend\nadd_Dantzig_Wolfe_decomposition(gap, dw_fct)Generalized Assignment problem can be solved using a Dantzig-Wolfe decomposition. Imagine we have a julia function that can solve efficiently the knapsack problem and returns the solution and the value of the solution(sol,value) = solveKnp(costs::Vector{Float64}, weights::Vector{Integer}, capacity::Integer)"
},

{
    "location": "callbacks.html#BlockDecomposition.getspid",
    "page": "Callbacks",
    "title": "BlockDecomposition.getspid",
    "category": "function",
    "text": "getspid(data::OracleSolverData)\n\nReturns the subproblem index for which the oracle has been assigned.\n\n\n\n"
},

{
    "location": "callbacks.html#BlockDecomposition.getcurcost",
    "page": "Callbacks",
    "title": "BlockDecomposition.getcurcost",
    "category": "function",
    "text": "getcurcost(x::JuMP.Variable)\n\nReturns the current cost of the varibale x.\n\n\n\n"
},

{
    "location": "callbacks.html#BlockDecomposition.getcurub",
    "page": "Callbacks",
    "title": "BlockDecomposition.getcurub",
    "category": "function",
    "text": "getcurub(x::JuMP.Variable)\n\nReturns the current ub of the varibale x.\n\n\n\n"
},

{
    "location": "callbacks.html#BlockDecomposition.getcurlb",
    "page": "Callbacks",
    "title": "BlockDecomposition.getcurlb",
    "category": "function",
    "text": "getcurlb(x::JuMP.Variable)\n\nReturns the current lb of the varibale x.\n\n\n\n"
},

{
    "location": "callbacks.html#BlockDecomposition.getcurdual",
    "page": "Callbacks",
    "title": "BlockDecomposition.getcurdual",
    "category": "function",
    "text": "getcurdual(x::JuMP.Constraint)\n\nReturns the current dual of the constraint x.\n\n\n\n"
},

{
    "location": "callbacks.html#JuMP.setsolutionvalue",
    "page": "Callbacks",
    "title": "JuMP.setsolutionvalue",
    "category": "function",
    "text": "setsolutionvalue!(data::OracleSolverData, x, val::Real)\n\nAssigns the value val to the variable x in the solution of the oracle solver\n\n\n\n"
},

{
    "location": "callbacks.html#BlockDecomposition.setsolutionbestobjval",
    "page": "Callbacks",
    "title": "BlockDecomposition.setsolutionbestobjval",
    "category": "function",
    "text": "setsolutionbestobjval!(data::OracleSolverData, objval::Real)\n\nAssigns the value value to the variable x in the solution of the oracle solver\n\n\n\n"
},

{
    "location": "callbacks.html#JuMP.addsolution",
    "page": "Callbacks",
    "title": "JuMP.addsolution",
    "category": "function",
    "text": "addsolution(data::OracleSolverData)\n\nIt ends the current solution and create a new solution in the oracle solver solution. Note that the previous solutions cannot be modified anymore.\n\n\n\n"
},

{
    "location": "callbacks.html#Write-the-oracle-solver-1",
    "page": "Callbacks",
    "title": "Write the oracle solver",
    "category": "section",
    "text": "We define an oracle that calls this function and solves each knapsack subproblem. ::function myKnapsackSolver(od::OracleSolverData)\n    machine = getspid(od)[0] # get the machine index\n    costs = [getcurcost(x[machine,j]) for j in Jobs] # get the current cost\n    (sol_x_m, value) = solveKnp(costs, Weight[m,:], Capacity[m]) # call the solver\n\n    # Building the oracle solution\n    for j in data.jobs\n        # add to oracle solution variables x[machine,j] with values sol_x_m[j]\n        setsolutionvalue(od, x[machine,j], sol_x_m[j])\n    end\n\n    # Set the objective value of the solution\n    setsolutionbestobjval(od, value)\nendgetspidgetcurcostgetcurubgetcurlbgetcurdualsetsolutionvaluesetsolutionbestobjvaladdsolution"
},

{
    "location": "callbacks.html#BlockDecomposition.addoracletosp!",
    "page": "Callbacks",
    "title": "BlockDecomposition.addoracletosp!",
    "category": "function",
    "text": "addoracletosp!(model::JuMP.Model, sp_id, sp_type::Symbol, f::Function)\n\nAttaches the oracle f to the subproblem of type sp_type which has the index spid. The argument spid must be a  Tuple or an  Integer.\n\n\n\n"
},

{
    "location": "callbacks.html#Attach-the-oracle-solver-1",
    "page": "Callbacks",
    "title": "Attach the oracle solver",
    "category": "section",
    "text": "Once the oracle solver function defined, we assign it to some subproblems using the following function.addoracletosp!In our example, for each subproblem with m, we assign the oracle myKnapsackSolver.for m in data.machines\n    addoracletosp!(gap, m, myKnapsackSolver)\nend"
},

{
    "location": "advanced.html#",
    "page": "Advanced",
    "title": "Advanced",
    "category": "page",
    "text": ""
},

{
    "location": "advanced.html#Advanced-features-1",
    "page": "Advanced",
    "title": "Advanced features",
    "category": "section",
    "text": ""
},

{
    "location": "advanced.html#BlockDecomposition.objectivevaluemagnitude",
    "page": "Advanced",
    "title": "BlockDecomposition.objectivevaluemagnitude",
    "category": "function",
    "text": "objectivevaluemagnitude(m::JuMP.Model, magnitude)\n\nSet the magnitude of the objective function of the model m to magnitude\n\n\n\n"
},

{
    "location": "advanced.html#BlockDecomposition.objectivevaluelowerbound",
    "page": "Advanced",
    "title": "BlockDecomposition.objectivevaluelowerbound",
    "category": "function",
    "text": "objectivevaluelowerbound(m::JuMP.Model, lb)\n\nSet the lower bound of the objective function of the model m to lb\n\n\n\n"
},

{
    "location": "advanced.html#BlockDecomposition.objectivevalueupperbound",
    "page": "Advanced",
    "title": "BlockDecomposition.objectivevalueupperbound",
    "category": "function",
    "text": "objectivevalueupperbound(m::JuMP.Model, ub)\n\nSet the upper bound of the objective function of the model m to ub\n\n\n\n"
},

{
    "location": "advanced.html#Objective-function-1",
    "page": "Advanced",
    "title": "Objective function",
    "category": "section",
    "text": "objectivevaluemagnitudeobjectivevaluelowerboundobjectivevalueupperbound"
},

{
    "location": "advanced.html#BlockDecomposition.branchingpriorityinmaster",
    "page": "Advanced",
    "title": "BlockDecomposition.branchingpriorityinmaster",
    "category": "function",
    "text": "branchingpriorityinmaster(x::JuMP.JuMPContainer, subproblem::Tuple{Symbol, Union{Tuple, Integer}}, priority)\n\nAssign to the variables x defined in the subproblem subproblem the priority value priority in master.\n\nbranchingpriorityinmaster(x, (:B_SP, 1), 2)\n\nThe variable x defined in the Benders subproblem with id 1 will have the branching priority value 2in the master.\n\n\n\n"
},

{
    "location": "advanced.html#BlockDecomposition.branchingpriorityinsubproblem",
    "page": "Advanced",
    "title": "BlockDecomposition.branchingpriorityinsubproblem",
    "category": "function",
    "text": "branchingpriorityinsubproblem(x::JuMP.JuMPContainer, subproblem::Tuple{Symbol, Union{Tuple, Integer}}, priority)\n\nAssign to the variables x defined in the subproblem subproblem the priority value priority in subproblems.\n\nbranchingpriorityinsubproblem(x, (:B_SP, 1), 2)\n\nThe variable x defined in the Benders subproblem with id 1 will have the branching priority value 2in subproblems.\n\n\n\n"
},

{
    "location": "advanced.html#BlockDecomposition.addbranching",
    "page": "Advanced",
    "title": "BlockDecomposition.addbranching",
    "category": "function",
    "text": "addbranching(model::JuMP.Model, rule::Symbol, varname::Symbol; args...)\n\ncreate a branching rule named rule on variable varname. Agruments are provided by the used and store in an array of pair. Arguments are checked by the solver.\n\n\n\n"
},

{
    "location": "advanced.html#Branching-priorities-1",
    "page": "Advanced",
    "title": "Branching priorities",
    "category": "section",
    "text": "branchingpriorityinmasterbranchingpriorityinsubproblemaddbranching"
},

{
    "location": "BlockSolverInterface.html#",
    "page": "BlockSolverInterface",
    "title": "BlockSolverInterface",
    "category": "page",
    "text": ""
},

{
    "location": "BlockSolverInterface.html#BlockSolverInterface-module-(dev)-1",
    "page": "BlockSolverInterface",
    "title": "BlockSolverInterface module (dev)",
    "category": "section",
    "text": "The content of this section is not useful to the user of BlockDecomposition. It has implementation details about the connection of BlockDecomposition with the underlying solver. It is aimed mainly at developpers who would like to contribute to BlockDecomposition."
},

{
    "location": "BlockSolverInterface.html#BlockDecomposition.BlockSolverInterface.set_constrs_decomposition!",
    "page": "BlockSolverInterface",
    "title": "BlockDecomposition.BlockSolverInterface.set_constrs_decomposition!",
    "category": "function",
    "text": "set_constrs_decomposition!(s::AbstractMathProgSolver, data::Array)\n\nsends to the solver s in which subproblems are the constraints. Each element of the data array is a Tuple containing (constr_name::Symbol, constr_id::Tuple, sp_type::Symbol, sp_id::Tuple). constr_name and constr_id are the name and the index of the constraint in the JuMP model. sp_type and sp_id are the type and the index of the subproblem to which the constraint is assigned.\n\n\n\n"
},

{
    "location": "BlockSolverInterface.html#BlockDecomposition.BlockSolverInterface.set_vars_decomposition!",
    "page": "BlockSolverInterface",
    "title": "BlockDecomposition.BlockSolverInterface.set_vars_decomposition!",
    "category": "function",
    "text": "set_vars_decomposition!(s::AbstractMathProgSolver, data::Array)\n\nsends to the solver s in which subproblems are the variables. Each element of the data array is a Tuple containing (var_name::Symbol, var_id::Tuple, sp_type::Symbol, sp_id::Tuple). var_name and var_id are the name and the index of the variable in the JuMP model. sp_type and sp_id are the type and the index of the subproblem to which the variable is assigned.\n\n\n\n"
},

{
    "location": "BlockSolverInterface.html#Decomposition-data-1",
    "page": "BlockSolverInterface",
    "title": "Decomposition data",
    "category": "section",
    "text": "BlockDecomposition.BlockSolverInterface.set_constrs_decomposition!BlockDecomposition.BlockSolverInterface.set_vars_decomposition!BlockDecomposition creates the decomposition list for both constraints and variables regardless of the type of decomposition used. Types of subproblem are ::DW_MASTER Dantzig-Wolfe master problem\n:B_MASTER Benders master problem\n:DW_SP Dantzig-Wolfe subproblem\n:B_SP Benders subproblem"
},

{
    "location": "BlockSolverInterface.html#BlockDecomposition.BlockSolverInterface.set_oracles!",
    "page": "BlockSolverInterface",
    "title": "BlockDecomposition.BlockSolverInterface.set_oracles!",
    "category": "function",
    "text": "set_oracles!(s::AbstractMathProgSolver, oracles::Array)\n\nsends to the solver s the list with subproblems and oracles functions. Each element of the oracles array is a Tuple containing (sp_id::Tuple, sp_type::Symbol, oracle::Function) with oracle, the function defined by the user.\n\n\n\n"
},

{
    "location": "BlockSolverInterface.html#BlockDecomposition.BlockSolverInterface.set_sp_mult!",
    "page": "BlockSolverInterface",
    "title": "BlockDecomposition.BlockSolverInterface.set_sp_mult!",
    "category": "function",
    "text": "set_sp_mult!(s::AbstractMathProgSolver, multiplicities::Array)\n\nsends to the solver s the multiplicity of each subproblem. Each element of the multiplicities array is a Tuple containing (sp_id::Tuple, sp_type::Symbol, mult_lb, mult_ub) with mult_lb the lower bound of the multiplicity and mult_ub the upper bound of the multiplicity.\n\n\n\n"
},

{
    "location": "BlockSolverInterface.html#BlockDecomposition.BlockSolverInterface.set_sp_prio!",
    "page": "BlockSolverInterface",
    "title": "BlockDecomposition.BlockSolverInterface.set_sp_prio!",
    "category": "function",
    "text": "set_sp_prio!(s::AbstractMathProgSolver, priorities::Array)\n\nsends to the solver s the list of subproblem priorities. Each element of the data array is a Tuple containing (sp_id::Tuple, sp_type::Symbol, sp_prio).\n\n\n\n"
},

{
    "location": "BlockSolverInterface.html#BlockDecomposition.BlockSolverInterface.set_var_branching_prio!",
    "page": "BlockSolverInterface",
    "title": "BlockDecomposition.BlockSolverInterface.set_var_branching_prio!",
    "category": "function",
    "text": "set_var_branching_prio!(s::AbstractMathProgSolver, priorities::Array)\n\nsends to the solver s the list of variables branching priorities. The number stored at the row i is the branching priority of the variable stored at the column i in the JuMP model.\n\n\n\n"
},

{
    "location": "BlockSolverInterface.html#BlockDecomposition.BlockSolverInterface.set_branching_rules!",
    "page": "BlockSolverInterface",
    "title": "BlockDecomposition.BlockSolverInterface.set_branching_rules!",
    "category": "function",
    "text": "set_branching_rules!(s::AbstractMathProgSolver, rules::Dict{Symbol, Any})\n\nsends to the solver s a Dictionary rules. The key is the name of the branching rule and the content is an array of branching instances. The array contains Tuple of variables name and parameters. Parameters are stored in an array of tuple (name_of_parameter, value_of_parameter).\n\nFor instance,\n\nrules = (:branching_rule_name => [(:x, [(:priority, 1)]), (:y, [(:priority, 1)])])\n\nNames of branching rules depend on solvers.\n\n\n\n"
},

{
    "location": "BlockSolverInterface.html#Additional-data-to-the-decomposition-1",
    "page": "BlockSolverInterface",
    "title": "Additional data to the decomposition",
    "category": "section",
    "text": "BlockDecomposition.BlockSolverInterface.set_oracles!BlockDecomposition.BlockSolverInterface.set_sp_mult!BlockDecomposition.BlockSolverInterface.set_sp_prio!BlockDecomposition.BlockSolverInterface.set_var_branching_prio!BlockDecomposition.BlockSolverInterface.set_branching_rules!"
},

{
    "location": "BlockSolverInterface.html#BlockDecomposition.BlockSolverInterface.set_objective_bounds_and_magnitude!",
    "page": "BlockSolverInterface",
    "title": "BlockDecomposition.BlockSolverInterface.set_objective_bounds_and_magnitude!",
    "category": "function",
    "text": "set_objective_bounds_and_magnitude!(s::AbstractMathProgSolver, magn, lb, ub)\n\nsends to the solver s the magnitude magn, the lower bound lb and the upper bound ub of the objective function.\n\n\n\n"
},

{
    "location": "BlockSolverInterface.html#Additional-data-to-the-model-1",
    "page": "BlockSolverInterface",
    "title": "Additional data to the model",
    "category": "section",
    "text": "BlockDecomposition.BlockSolverInterface.set_objective_bounds_and_magnitude!Send to the solver s the magnitude magn, the lower bound lb and the upper bound ub of the objective function."
},

{
    "location": "BlockSolverInterface.html#BlockDecomposition.BlockSolverInterface.getcurrentcost",
    "page": "BlockSolverInterface",
    "title": "BlockDecomposition.BlockSolverInterface.getcurrentcost",
    "category": "function",
    "text": "getcurrentcost(m::AbstractMathProgModel, vcol::Integer)\n\nreturns the current cost of the vcol ^th variable.\n\n\n\n"
},

{
    "location": "BlockSolverInterface.html#BlockDecomposition.BlockSolverInterface.getdisaggregatedvalueofvariable",
    "page": "BlockSolverInterface",
    "title": "BlockDecomposition.BlockSolverInterface.getdisaggregatedvalueofvariable",
    "category": "function",
    "text": "getdisaggregatedvalueofvariable(m::AbstractMathProgModel, vcol::Integer)\n\nreturns the disaggregated value of the vcol ^th variable.\n\n\n\n"
},

{
    "location": "BlockSolverInterface.html#Costs-and-solutions-1",
    "page": "BlockSolverInterface",
    "title": "Costs and solutions",
    "category": "section",
    "text": "BlockDecomposition.BlockSolverInterface.getcurrentcostBlockDecomposition.BlockSolverInterface.getdisaggregatedvalueofvariable"
},

{
    "location": "BlockSolverInterface.html#BlockDecomposition.BlockSolverInterface.set_oraclesolution_solution",
    "page": "BlockSolverInterface",
    "title": "BlockDecomposition.BlockSolverInterface.set_oraclesolution_solution",
    "category": "function",
    "text": "set_oraclesolution_solution(o::OracleSolverData, x::JuMP.Variable, v::Real)\n\nSet the value of the variable x to v in the oracle solver solution stored in the OracleSolverData object o.\n\n\n\n"
},

{
    "location": "BlockSolverInterface.html#BlockDecomposition.BlockSolverInterface.set_oraclesolution_objval",
    "page": "BlockSolverInterface",
    "title": "BlockDecomposition.BlockSolverInterface.set_oraclesolution_objval",
    "category": "function",
    "text": "set_oraclesolution_objval(o::OracleSolverData, v::Real)\n\nSets the objective value stored in the OracleSolverData object o to v.\n\n\n\n"
},

{
    "location": "BlockSolverInterface.html#BlockDecomposition.BlockSolverInterface.set_oraclesolution_newsolution",
    "page": "BlockSolverInterface",
    "title": "BlockDecomposition.BlockSolverInterface.set_oraclesolution_newsolution",
    "category": "function",
    "text": "set_oraclesolution_newsolution(o::OracleSolverData)\n\ncreates a new solution in the oracle solver solution. It is usefull, if the user wants to return several solutions.\n\n\n\n"
},

{
    "location": "BlockSolverInterface.html#BlockDecomposition.BlockSolverInterface.get_oracle_phaseofstageapproach",
    "page": "BlockSolverInterface",
    "title": "BlockDecomposition.BlockSolverInterface.get_oracle_phaseofstageapproach",
    "category": "function",
    "text": "get_oracle_phaseofstageapproach(o::OracleSolverData)\n\nReturns the phase of stage approach.\n\n\n\n"
},

{
    "location": "BlockSolverInterface.html#Oracle-solver-1",
    "page": "BlockSolverInterface",
    "title": "Oracle solver",
    "category": "section",
    "text": "BlockDecomposition.BlockSolverInterface.set_oraclesolution_solutionBlockDecomposition.BlockSolverInterface.set_oraclesolution_objvalBlockDecomposition.BlockSolverInterface.set_oraclesolution_newsolutionBlockDecomposition.BlockSolverInterface.get_oracle_phaseofstageapproach"
},

]}

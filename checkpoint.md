---
title: Project Checkpoint
permalink: checkpoint.html
---

# Project Checkpoint

### Schedule

- 4/24 ~ 4/27
  - Conduct OpenMP experiments on Matrix and Vector manipulation functions to
    decide at which size we should enable OpenMP parallelism.
  - Implemented OpenMP version functions for Matrix and Vector manipulations.
- 4/28 ~ 4/30
  - Optimize memory usage for Matrix and Vector manipulations.
  - Conduct CUDA experiments for Matrix and Vector manipulations.
- 5/1 ~ 5/4
  - Implement CUDA version Matrix and Vector manipulations.
  - Survey possible parallelization in the Krylov process.
- 5/5 ~ 5/7
  - Buffer time; parallelize the Krylov process.
- 5/8 ~ 5/10
  - Conduct benchmarks and comparisons.
- 5/10 ~ 5/12
  - Buffer time; organize results and complete the project page.
  
  
### Progress (What we have done so far)

We implemented a sequential version of the GMRES algorithm, and we also built
our own sparse matrix structure (Compressed Sparse Row storage). We compared the
dense and sparse implementations. Also, we experimented at what size the benefit
of using OpenMP will overcome the overhead. We also optimized our data structure
to make it more easily to be parallelized.

<!-- The checkpoint exists is to give you a deadline approximately halfway through -->
<!-- the project. The following are suggestions for information to include in your -->
<!-- checkpoint write-up. Your goal in the writeup is to assure the course staff (and -->
<!-- yourself) that your project is proceeding as you said it would in your proposal. -->
<!-- If it is not, your checkpoint writeup should emphasize what has been causing you -->
<!-- problems, and provide an adjusted schedule and adjusted goals. As projects -->
<!-- differ, not all items in the list below are relevant to all projects. -->

<!-- - Make sure your project schedule on your main project page is up to date with -->
<!--   work completed so far, and well as with a revised plan of work for the coming -->
<!--   weeks. As by this time you should have a good understanding of what is -->
<!--   required to complete your project, I want to see a very detailed schedule for -->
<!--   the coming weeks. I suggest breaking time down into half-week increments. Each -->
<!--   increment should have at least one task, and for each task put a person's name -->
<!--   on it. -->
<!-- - One to two paragraphs, summarize the work that you have completed so far. -->
<!--   (This should be easy if you have been maintaining this information on your -->
<!--   project page.) -->
<!-- - Describe how you are doing with respect to the goals and deliverables stated -->
<!--   in your proposal. Do you still believe you will be able to produce all your -->
<!--   deliverables? If not, why? What about the "nice to haves"? In your checkpoint -->
<!--   writeup we want a new list of goals that you plan to hit for the Parallelism -->
<!--   competition. -->
<!-- - What do you plan to show at the parallelism competition? Will it be a demo? -->
<!--   Will it be a graph? -->
<!-- - Do you have preliminary results at this time? If so, it would be great to -->
<!--   included them in your checkpoint write-up. -->
<!-- - List the issues that concern you the most. Are there any remaining unknowns -->
<!--   (things you simply don't know how to solve, or resource you don't know how to -->
<!--   get) or is it just a matter of coding and doing the work? If you do not wish -->
<!--   to put this information on a public web site you are welcome to email the -->
<!--   staff directly. -->
<!-- - [Optionally] schedule a meeting with the course staff to discuss progress -->


### Deliverable revisit

We list four deliverables in our project proposal: sequential version, OpenMP 
version, GPU version, and CPU-GPU mixed version using OpenMPI. We planned to 
complete GPU version before the project checkpoint but we are not able to meet 
the milestones. We finished the sequential version. Currently, we are working 
on the OpenMP version, which we may finish in two to three days after checkpoint; 
that is, we have a delay about a week. We may not complete all planned deliverable. 
Because the GPU version requires a re-design of data structure, we may not have 
enough time for CPU-GPU mixed version across multiple nodes. We decide to move 
this part into nice to have features. By the deadline of final project, we should 
be able to complete GPU version. Our new milestones will be shifted by a week. 
We anticipa e to complete a workable GPU version before May 3. 

### Issues

We encounter a problem related to solver precision while implementing the sequential 
version. Some large matrices do not converge so that we can not use them to compare 
with the OpenMPI version. We spent a lot of time to survey this convergence issue. 
We finally realized that GMRES algorithm requires proper precondition to reduce 
condition numbers. This reduction is highly matrix dependent. We found that different 
domain may have different precondition matrices. Therefore, we try to find some matrix 
patterns that does not rely on precondition and is likely to converge as long as the 
number of iterations is enough. We found this kind of matrices on [Link]. Now, we can 
focus on how to parallelize algorithm instead of being bothered by algorithm itself. 
As for environment configuration, we think we don’t have problem so far.

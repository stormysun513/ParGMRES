Parallel GMRES: A parallel linear solver.
============

## TITLE

Parallel GMRES: A parallel linear solver.

## TEAM MEMBER

- Yu-Lun Tsai (yulunt@andrew)
- Chih-Wei Chang (cchang3@andrew)

## SUMMARY

We propose to implement a parallelized GMRES (Generalized Minimal RESidual
method), which is a commonly used solver for large sparse linear system, on
GPU-CPU heterogeneous platform.

## BACKGROUND

We try to build a parallel version implementation of GMRES. GMRES is one of the
mostly used algorithm in solving large sparse linear system. Solving a linear
system efficiently can be crucial in that many computate-intensive applications,
such as machine learning or scientific computation, involve solving different
kind of linear systems. GMRES is by itself an iterative method, in which the
main flow of the algorithm is inherently sequential. However, the Arnoldi
process, the most time-consuming part in GMRES, is consisted of generating basis
and orthogonalization, which is potentially paralleizable. Hence, we belive a
well-designed Parallel GMRES should have the chance to out-perform the
sequential GMRES by a significant scale.

<!-- If your project involves accelerating a compute-intensive application, describe -->
<!-- the application or piece of the application you are going to implement in more -->
<!-- detail. This description need only be a few paragraphs. It might be helpful to -->
<!-- include a block diagram or pseudocode of the basic idea. An important detail is -->
<!-- what aspects of the problem might benefit from parallelism? and why? -->


## THE CHALLENGE

Describe why the problem is challenging. What aspects of the problem might make
it difficult to parallelize? In other words, what to you hope to learn by doing
the project?

- Describe the workload: what are the dependencies, what are its memory access
  characteristics? (is there locality? is there a high communication to
  computation ratio?), is there divergent execution?
- Describe constraints: What are the properties of the system that make mapping
  the workload to it challenging?

## RESOURCES
Describe the resources (type of computers, starter code, etc.) you will use.
What code base will you start from? Are you starting from scratch or using an
existing piece of code? Is there a book or paper that you are using as a
reference (if so, provide a citation)? Are there any other resources you need,
but haven't figured out how to obtain yet? Could you benefit from access to any
special machines?

## GOALS AND DELIVERABLES

This is by far the most important section of the proposal:

- Separate your goals into what you PLAN TO ACHIEVE (what you believe you must
  get done to have a successful project and get the grade you expect) and an
  extra goal or two that you HOPE TO ACHIEVE if the project goes really well and
  you get ahead of schedule. It may not be possible to state precise performance
  goals at this time, but we encourage you be as precise as possible. If you do
  state a goal, give some justification of why you think you can achieve it.
  (e.g., I hope to speed up my starter code 10x, because if I did it would run
  in real-time.)
- If applicable, describe the demo you plan to show at the parallelism
  computation (will it be an interactive demo? will you show an output of the
  program that is really neat? will you show speedup graphs?). Specifically,
  what will you show us that will demonstrate you did a good job?
- If your project is an analysis project, what are you hoping to learn about the
  workload or system being studied? What question(s) do you plan to answer in
  your analysis?
- Systems project proposals should describe what the system will be capable of
  and what performance is hoped to be achieved.
- IN GENERAL: Imagine that I didn't give you a grading script on assignments 2,
  3, or 4. Imagine you did the entire assignment, made it as fast as you could,
  and then turned it in. You wouldn't have any idea if you'd done a good job!!!
  That's the situation you are in for the final project. And that's the
  situation I'm in when grading your final project. As part of your project
  plan, and ONE OF THE FIRST THINGS YOU SHOULD DO WHEN YOU GET STARTED WORKING
  is implement the test harnesses and/or baseline "reference" implementations
  for your project. Then, for the rest of your project you always have the
  ability to run your optimized code and obtain a comparison.

## PLATFORM CHOICE
Describe why the platform (computer and/or language) you have chosen is a good
one for your needs. Why does it make sense to use this parallel system for the
workload you have chosen?

## SCHEDULE
Produce a schedule for your project. Your schedule should have at least one item
to do per week. List what you plan to get done each week from now until the
parallelism competition in order to meet your project goals. Keep in mind that
due to other classes, you'll have more time to work some weeks than others (work
that into the schedule). You will need to re-evaluate your progress at the end
of each week and update this schedule accordingly. Note the intermediate
checkpoint deadline is April 25th. In your schedule we encourage you to be
precise as precise as possible. It's often helpful to work backward in time from
your deliverables and goals, writing down all the little things you'll need to
do (establish the dependencies!).

## Reference

[1] I am reference 1

<!-- ## Welcome to GitHub Pages -->

<!-- You can use the [editor on GitHub](https://github.com/stormysun513/pcap-final/edit/gh-pages/README.md) to maintain and preview the content for your website in Markdown files. -->

<!-- Whenever you commit to this repository, GitHub Pages will run [Jekyll](https://jekyllrb.com/) to rebuild the pages in your site, from the content in your Markdown files. -->

<!-- ### Markdown -->

<!-- Markdown is a lightweight and easy-to-use syntax for styling your writing. It includes conventions for -->

<!-- ```markdown -->
<!-- Syntax highlighted code block -->

<!-- # Header 1 -->
<!-- ## Header 2 -->
<!-- ### Header 3 -->

<!-- - Bulleted -->
<!-- - List -->

<!-- 1. Numbered -->
<!-- 2. List -->

<!-- **Bold** and _Italic_ and `Code` text -->

<!-- [Link](url) and ![Image](src) -->
<!-- ``` -->

<!-- For more details see [GitHub Flavored Markdown](https://guides.github.com/features/mastering-markdown/). -->

<!-- ### Jekyll Themes -->

<!-- Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/stormysun513/pcap-final/settings). The name of this theme is saved in the Jekyll `_config.yml` configuration file. -->

<!-- ### Support or Contact -->

<!-- Having trouble with Pages? Check out our [documentation](https://help.github.com/categories/github-pages-basics/) or [contact support](https://github.com/contact) and we’ll help you sort it out. -->

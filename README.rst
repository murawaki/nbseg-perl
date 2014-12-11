Obsolete Perl code for non-parametric Bayesian text segmentation
================================================================================

This repository contains Perl code I no longer use, including

* a group of Dirichlet/Pitman-Yor processes,
* a character-based bigram word model, and
* unigram/bigram word models with token-based, block and type-based sampling.


Requirements
------------------------------------------------------------

The following CPAN modules are required:

* Math::GSL
* MAth::Cephes
* Regexp::Assemble
* Carp::Assert



Run a sample script
------------------------------------------------------------

% perl -Ilib scripts/sample-token.pl --seed=1 --type=Dirichlet --input=samples/alice.unseg --iter=100 --nested --debug --randInit=0.1

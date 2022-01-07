# Digital Communications Project-2021

## Convolutional Coding and the Viterbi Algorithm

------------------------------------------------------------------------------------------------------------

**Task graphs/code are available to view by checking out the corresponding task branch**

ex. `git checkout task-1`

Alternatively, if one does not have git installed, the respective `main.m` files are stored in `/matlab/main-archives/main-task-n`.

------------------------------------------------------------------------------------------------------------

## Description of scripts & classes:

#### main.m

Main script, used to launch the simulation chain and output the figures.

#### ConvEncoder.m

Provides the description and functionalities of all the different convolutional encoders.

#### DecoderType.m

Defines the two types [HARD/SOFT] of Viterbi decoding algorithms and their implementation.

#### ViterbiDecoder.m

Modulare handle for storing a decoder object with its given decoder type and other properties

#### Trellis.m

Definitions and implementations of all the different encoder trellis

#### SymbolMapper.m

Contains all implementations of the constellations, as well as helper methods like the AWGN channel simulation.

#### LoadingBar.m

Simple class to display loading bar when running.

------------------------------------------------------------------------------------------------------------

#### Project Overview

The goal of this project is to provide a deeper understanding of how channel coding can be used to improve the performance of a digital communication system. The main task is to implement, evaluate, and compare several options for different encoders and modulation formats. In order to successfully complete this task, we have to program certain algorithms (like the Viterbi algorithm), as well as answer specific questions with the help of concepts introduced in the lectures (like channel capacity).The deliverables of the project are split into two parts.

#### Learning Outcomes 

After completing this project, we should be able to

- successfully set up and run a simulation environment efficiently in MATLAB,
- understand why and how MATLAB can be used to simulate the performance of a system that eventually transmits and receives analog waveforms
- use the concept of channel capacity as a benchmarking tool 
- encode a bit stream with a convolutional encoder and understand the relationship between the encoder and a trellis diagram
- implement the Viterbi algorithm and understand how it works
- explain and quantify the performance difference between hard and soft decoding
- analyze fundamental trade-offs between power, bit rate, bandwidth, and bit-error rate in a coded digital communication system and in particular make fair comparisons between uncoded and coded systems.

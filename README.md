# Julia for Bioinformaticians
Code for special course in Julia aug-sep 2020

## Project (QT Clustering)
### [Report: GitHub Wiki](https://github.com/FSGade/julia_course/wiki/Report:-Quality-Threshold-Clustering)

>QT (Quality Threshold) Clustering is an algorithm that groups multi-dimensional vectors into high quality clusters. Quality is ensured by finding large cluster whose diameter does not exceed a given user-defined diameter threshold.
>This method prevents dissimilar vectors from being forced under the same cluster and ensures that only good quality clusters will be formed.

```
usage: julia qt.jl inputfile threshold

arguments:
    inputfile	Input file name (tab-separated list of floats w/ optional name
		in first column)
    threshold	Quality threshold of cluster diameter
```

## Exercises
The exercise script files can be found in exercises/ex*n*.jl

| n | Exercise name | Link |
| --- | ------------- | ----------- |
| 1 | ~Python~ Julia Basics |  [[1]](https://teaching.healthtech.dtu.dk/22110/index.php/Python_Basics)|
| 2 | ~Python~ Julia Simple file reading |  [[2]](https://teaching.healthtech.dtu.dk/22110/index.php/Python_Simple_file_reading)|
| 3 | ~Python~ Julia Input-Output |  [[3]](https://teaching.healthtech.dtu.dk/22110/index.php/Python_Input-Output)|
| 4 | Exceptions and Bug Handling |  [[4]](https://teaching.healthtech.dtu.dk/22110/index.php/Exceptions_and_Bug_Handling)|
| 5 | Lists/Sequences |  [[5]](https://teaching.healthtech.dtu.dk/22110/index.php/Lists/Sequences)|
| 6 | Pattern Matching and Regular Expressions |  [[6]](https://teaching.healthtech.dtu.dk/22110/index.php/Pattern_Matching_and_Regular_Expressions)|
| 7 | Sets and Dictionaries |  [[7]](https://teaching.healthtech.dtu.dk/22110/index.php/Sets_and_Dictionaries)|
| 8 | ~Python~ Julia and Advanced Data Structures |  [[8]](https://teaching.healthtech.dtu.dk/22110/index.php/Python_and_Advanced_Data_Structures)|
| 9 | Comprehension and Generators |  [[9]](https://teaching.healthtech.dtu.dk/22110/index.php/Comprehension_and_Generators)|
| 10 | Useful Functions and Methods |  [[10]](https://teaching.healthtech.dtu.dk/22110/index.php/Useful_Functions_and_Methods)|

**NB:** This does not include the [Python exercise on functions](https://teaching.healthtech.dtu.dk/22110/index.php/Python_Functions), as they are already implemented in the other exercises per the Julia idiom.

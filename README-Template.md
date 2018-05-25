# Fit-Seq

One Paragraph of project description goes here

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

What things you need to install the software and how to install them

```
Give examples
```

### Installing

A step by step series of examples that tell you have to get a development env running

Say what the step will be

```
Give the example
```

And repeat

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo

## Running the tests

Explain how to run the automated tests for this system

### Pooled growth simulation
A numerical simulation framework is developed to simulate the competitive pooled cell growth of a population of genotypes with different fitnesses, and includes all sources of experimental noise, such as growth noise, sampling during bottlenecks, DNA extraction, PCR, and sampling on the sequencer. The simulation is performed in MATLAB using “evolution_0mut_complex.m”.


#### Inputs
-- lineage: number of genotypes of the population

-- t_evo: total number of growth generations

-- cell_num_ini: a vector of initial cell number of each genotype at the 0-th generation, size = lineage * 1

-- x_ini: a vector of the fitness of each genotype, size = lineage * 1

-- deltat: number of generations between successive cell transfers

-- read_depth_average: average number of reads per genotype per sequencing time point

-- noise_option: a vector of options of whether five types of noise, cell growth, cell transfer at the bottleneck, DNA extraction, PCR, and sequencing (in that order in the vector) are simulated, size = 1 * 5, 1 means the noise is included and 0 means the noise is not included.

-- 'format': optional, file format of the output file, 'csv'(default) or 'mat'

#### Outputs
-- file_name: name of the file generated, the standard file generated 'data_evo_simu_0mut_complex_********-*********.mat' when 'format' is set to be 'mat', and 'data_evo_simu_0mut_complex_********-*********.csv' when 'format' is set to be 'csv'

```
lineage = 1e4;
t_evo = 24;
cell_num_ini = 1e2*ones(lineage,1);  
x_ini = random('Normal',0,0.15, [lineage,1]);  % Gaussian with mean = 0 and standard devistion = 0.15
read_depth_average = 100;  
deltat = 8;  
noise_option = [1,1,1,1,1];  % All five types of noise are simulated

[ file_name ] = evolution_0mut_simplified(lineage, t_evo, cell_num_ini, x_ini, noise_option);
```

#### Complex version

### Fitness estimation using Fit-Seq

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone who's code was used
* Inspiration
* etc


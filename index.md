
# A New Study of Segregation and Inequality

Demographic data is becoming more detailed and more accessible. The new availability of high-quality data has begun to make its mark in quantitative segregation studies. [Classical](https://academic.oup.com/sf/article/67/2/281/2231999/The-Dimensions-of-Residential-Segregation) [studies](http://journals.sagepub.com/doi/abs/10.1111/1467-9531.00110) in the field focused on a single question: 

> How (**to what degree**) segregated is this city? 

In the last few years, however, researchers have begun to ask a more complex question: 

> How (**in what way**) is this city segregated? 

This question has multiple layers: 

1. On what *scale* is the city segregated? Are there massive, homogeneous regions of race or income group? Or are there smaller pockets of demographic variation? While previous studies ([example](http://www.sciencedirect.com/science/article/pii/S0022053110001353)) have addressed the idea of hierarchical "levels" in segregation, a fully spatial treatment of the scale of segregation requires new methodology. 
2. How can we *visualize* segregation? Visualization is straightforward in the case of relatively few racial groups, as beautifully demonstrated by [this fabulous example](https://demographics.virginia.edu/DotMap/index.html). But when there are many groups, or when the groups have an intrinsic ordering, visualization may be more complex. 
3. How can we *find* the segregated structure of cities? This may be viewed as a problem of machine-learning or statistical inference, and has application for visualization and dimension-reduction purposes.  

# Introducing `compx`

Modern questions require modern tools. `compx` is an `R` package for analyzing segregation and inequality using contemporary mathematics and machine learning. The theoretical foundation of `compx` is [*information geometry*](https://en.wikipedia.org/wiki/Information_geometry), a relatively young branch of mathematics that uses the classical tools of differential geometry to learn about probabilistic and statistical models. Our core approach is to treat a map with demographics as a parameterized statistical model, and then use elementary notions of information geometry to learn about it. This approach enables the following two main forms of analysis. 

## Studying Spatial Complexity

How complex are the patterns of segregation in different cities? Where are the spatial "boundaries" between demographic groups? `compx` provides mathematical methods for approaching these questions through the unifying idea of the metric tensor. Learn more about computing, analyzing, and visualizing the metric tensor [here](https://philchodrow.github.io/compx/vignette_metric.html). 

## Spatial Network Analysis

What are the characteristic scales of segregation in a city? How can I visualize these scales? What are the natural spatial clusters from which the city is assembled? `compx` provides tools to approach these questions through network analysis. The core idea is to compute the [*geodesic distance*](https://en.wikipedia.org/wiki/Geodesic) between spatial locales, and use this distance as a graph weighting. Various forms of network analysis--such as spectral analysis and community-detection-- can then be used to characterize the spatial structure of segregation. Learn more about constructing, analyzing, and visualizing these networks [here](https://philchodrow.github.io/compx/clustering.html). 

# Questions, Help, Feedback

I'd love to hear from you! Please reach out through GitHub, or my [website](https://philchodrow.github.io/). 

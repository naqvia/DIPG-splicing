

# A plotting R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800


# --------------------------------------------------------------------------
# If you don't have the ggplot2 package installed, uncomment the following line:
# install.packages( "ggplot2" );
library( ggplot2 );
# --------------------------------------------------------------------------
# If you don't have the scales package installed, uncomment the following line:
# install.packages( "scales" );
library( scales );


# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","log10_p_value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0005925","focal adhesion", 2.125, 6.006,-3.170, 2.595,-3.3820,0.856,0.000),
c("GO:0005938","cell cortex", 1.287, 4.394, 4.229, 2.378,-8.3862,0.750,0.000),
c("GO:0016020","membrane",50.676, 4.907,-0.577, 3.972,-4.1244,0.972,0.000),
c("GO:0030054","cell junction", 6.490, 1.223,-1.393, 3.080,-3.4089,0.947,0.000),
c("GO:0044456","synapse part", 3.645,-1.076, 7.378, 2.829,-5.9066,0.921,0.000),
c("GO:0099080","supramolecular complex", 4.889,-0.588, 2.373, 2.957,-3.7959,0.946,0.000),
c("GO:0097458","neuron part", 7.069,-6.415,-5.224, 3.117,-7.0846,0.943,0.003),
c("GO:0044463","cell projection part", 5.425, 7.083, 0.513, 3.002,-5.9872,0.942,0.004),
c("GO:0042995","cell projection", 9.968,-3.794, 6.975, 3.266,-5.3675,0.944,0.004),
c("GO:0035371","microtubule plus-end", 0.103,-6.353, 2.195, 1.301,-4.3134,0.675,0.012),
c("GO:0005737","cytoplasm",60.330,-4.664,-4.071, 4.048,-5.5317,0.944,0.027),
c("GO:0099568","cytoplasmic region", 1.601, 3.465,-4.962, 2.473,-6.9666,0.898,0.083),
c("GO:0030139","endocytic vesicle", 1.476,-0.799,-6.584, 2.438,-4.3206,0.837,0.084),
c("GO:0005829","cytosol",25.879, 2.421,-6.311, 3.680,-4.4750,0.891,0.128),
c("GO:0031256","leading edge membrane", 0.741, 2.460, 5.607, 2.140,-3.4401,0.842,0.187),
c("GO:0005856","cytoskeleton",11.022,-5.462, 1.988, 3.309,-5.3261,0.765,0.267),
c("GO:0012506","vesicle membrane", 4.008,-1.933,-5.313, 2.870,-3.1824,0.824,0.307),
c("GO:0005886","plasma membrane",28.794, 3.127, 3.545, 3.726,-4.7825,0.836,0.309),
c("GO:0098590","plasma membrane region", 5.165, 2.644, 4.638, 2.980,-3.7423,0.824,0.426),
c("GO:0072686","mitotic spindle", 0.373,-6.691, 1.233, 1.845,-4.2125,0.726,0.427),
c("GO:0098835","presynaptic endocytic zone membrane", 0.005, 0.528, 6.598, 0.301,-3.9066,0.849,0.455),
c("GO:0030659","cytoplasmic vesicle membrane", 3.921,-1.132,-5.354, 2.861,-3.3420,0.765,0.508),
c("GO:0099513","polymeric cytoskeletal fiber", 3.678,-6.186, 0.851, 2.833,-3.6364,0.586,0.566),
c("GO:1990752","microtubule end", 0.130,-6.511, 1.813, 1.398,-3.7258,0.670,0.579),
c("GO:0030130","clathrin coat of trans-Golgi network vesicle", 0.065,-1.521,-6.343, 1.114,-3.2291,0.810,0.584),
c("GO:0045180","basal cortex", 0.027, 5.086, 4.416, 0.778,-4.4989,0.788,0.647),
c("GO:0030055","cell-substrate junction", 2.169, 5.823,-3.441, 2.604,-3.2916,0.872,0.666),
c("GO:0044430","cytoskeletal part", 8.345,-5.929, 0.461, 3.189,-7.9355,0.683,0.678),
c("GO:0005815","microtubule organizing center", 3.553,-5.977, 1.090, 2.818,-3.4855,0.672,0.686));

one.data <- data.frame(revigo.data);
names(one.data) <- revigo.names;
one.data <- one.data [(one.data$plot_X != "null" & one.data$plot_Y != "null"), ];
one.data$plot_X <- as.numeric( as.character(one.data$plot_X) );
one.data$plot_Y <- as.numeric( as.character(one.data$plot_Y) );
one.data$plot_size <- as.numeric( as.character(one.data$plot_size) );
one.data$log10_p_value <- as.numeric( as.character(one.data$log10_p_value) );
one.data$frequency <- as.numeric( as.character(one.data$frequency) );
one.data$uniqueness <- as.numeric( as.character(one.data$uniqueness) );
one.data$dispensability <- as.numeric( as.character(one.data$dispensability) );
#head(one.data);


# --------------------------------------------------------------------------
# Names of the axes, sizes of the numbers and letters, names of the columns,
# etc. can be changed below

p1 <- ggplot( data = one.data );
p1 <- p1 + geom_point( aes( plot_X, plot_Y, colour = log10_p_value, size = plot_size), alpha = I(0.6) ) + scale_size_area();
p1 <- p1 + scale_colour_gradientn( colours = c("blue", "green", "yellow", "red"), limits = c( min(one.data$log10_p_value), 0) );
p1 <- p1 + geom_point( aes(plot_X, plot_Y, size = plot_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) )) + scale_size_area();
p1 <- p1 + scale_size( range=c(5, 30)) + theme_bw(); # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
ex <- one.data [ one.data$dispensability < 0.15, ]; 
p1 <- p1 + geom_text( data = ex, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 3 );
p1 <- p1 + labs (y = "semantic space x", x = "semantic space y");
p1 <- p1 + theme(legend.key = element_blank()) ;
one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y);
p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/10,max(one.data$plot_X)+one.x_range/10);
p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/10,max(one.data$plot_Y)+one.y_range/10);



# --------------------------------------------------------------------------
# Output the plot to screen

p1;

# Uncomment the line below to also save the plot to a file.
# The file type depends on the extension (default=pdf).

# ggsave("C:/Users/path_to_your_file/revigo-plot.pdf");

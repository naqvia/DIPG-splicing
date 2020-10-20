

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
revigo.data <- rbind(c("GO:0001817","regulation of cytokine production", 0.108,-3.977, 4.967, 4.141,-4.9318,0.726,0.000),
c("GO:0002376","immune system process", 0.600, 3.093, 5.647, 4.886,-9.0799,0.953,0.000),
c("GO:0006952","defense response", 0.568,-4.764,-3.861, 4.863,-8.5157,0.631,0.000),
c("GO:0032940","secretion by cell", 0.763, 4.609,-2.963, 4.991,-4.2426,0.727,0.070),
c("GO:0001775","cell activation", 0.171, 5.128, 1.428, 4.341,-3.4895,0.858,0.166),
c("GO:0043900","regulation of multi-organism process", 0.079,-2.233, 4.605, 4.008,-4.1772,0.792,0.189),
c("GO:0008284","positive regulation of cell proliferation", 0.151,-5.460, 3.777, 4.288,-4.1152,0.629,0.199),
c("GO:0019221","cytokine-mediated signaling pathway", 0.093,-5.571,-1.324, 4.078,-4.2749,0.576,0.304),
c("GO:0032101","regulation of response to external stimulus", 0.160,-4.432,-0.023, 4.313,-4.9355,0.515,0.318),
c("GO:0006955","immune response", 0.337,-2.688,-4.014, 4.635,-7.1568,0.299,0.340),
c("GO:0007186","G-protein coupled receptor signaling pathway", 0.882,-4.499,-0.910, 5.054,-4.0496,0.552,0.433),
c("GO:0006954","inflammatory response", 0.110,-4.731,-4.595, 4.151,-6.6289,0.605,0.491),
c("GO:0046903","secretion", 0.810, 4.390,-4.285, 5.017,-3.8539,0.808,0.557),
c("GO:0032103","positive regulation of response to external stimulus", 0.055,-5.085, 0.570, 3.845,-4.5834,0.494,0.637),
c("GO:0097305","response to alcohol", 0.055,-6.109,-3.365, 3.850,-3.0367,0.680,0.645),
c("GO:0043299","leukocyte degranulation", 0.014, 3.129,-4.419, 3.252,-3.2840,0.363,0.698));

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

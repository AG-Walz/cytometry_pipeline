# Create a frame to hold multiple plots
par(mfrow = c(3, 7))

# Plot sequentially using a loop
plots <- lapply(1:21, function(i) {
  flowViz::plot(flow_frames[[i]], "CD69")
})
flowViz::plot(flow_frames[[1]])

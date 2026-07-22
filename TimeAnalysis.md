On the standard settings, with curvature = $800$, there are $40062$ target points.

The full run time was $858$ seconds

The time to calculate the quad volume potential contributions was $798$ seconds.

Of this time, the time spent calculating F matrices was $1.8$ seconds.

The time spent calculating L matrices was $60.2$ seconds.

The time separating meshes was $666.5$ seconds (!!!)

Now, a table describing the number and percent of cache hits:

| Total cache hits | F cache hit number | F cache hit % | L cache hit number | L cache hit % |
| ---------------- | ------------------ | ------------- | ------------------ | ------------- |
| $185680$         | $183572$           | $98.8$%       | $130906$           | $70.5$%       |

The time spent calculating L matrices can probably be cut down.

The main time gain will undoubtedly be fixing the separateMesh function.
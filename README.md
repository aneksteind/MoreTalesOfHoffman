Code written in support of testing various bounds in [More Tales of Hoffman: bounds for the vector chromatic number of a graph](https://arxiv.org/abs/1812.02613)

## run instructions
```
python conj.py test --help
Usage: conj.py test [OPTIONS]

Options:
  --max INTEGER                   the maximum number of vertices  [required]
  --min INTEGER                   the minimum number of vertices
  --graph-source [named|circulant]
                                  named graphs from the mathematica database
                                  or generated circulant graphs  [required]
  --out TEXT                      the .csv file to write results to
                                  [required]
  --verbose / --quiet             print the graph names during execution
  --help                          Show this message and exit.
```
  
## caveats
Running this code with the intention of using `named` graphs requires a directory of the name `GraphData`, with DOT files of graphs contained within. More about DOT files can be found [here](https://en.wikipedia.org/wiki/DOT_%28graph_description_language%29).

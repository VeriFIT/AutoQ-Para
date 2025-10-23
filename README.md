# AutoQ-Para 
A tool for parametrized verification of quantum circuits using the Synchronized Weighted Tree Automata (SWTAs)

## Building and running
In order to streamline the build process, we provide a container file
allowing to build and execute the tool on a predefined set of verification
problems.

```bash
podman build -t swta-impl .
podman run --rm swta-impl
```

To see the list of dependencies for modifying the tool, see `Dockerfile`.

## Publications
- POPL'26 (Link: TBD)

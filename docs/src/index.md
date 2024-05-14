# AlgebraicABMs.jl

```@meta
CurrentModule = AlgebraicABMs
```

`AlgebraicABMs.jl` is a Julia library for creating agent-based models. We ultimate want to provide capabilities on par with software like AnyLogic, NetLogo, and Agents.jl - all while offering a mostly declarative interface such that: 
- The model's logic is transparent (rather than hidden away in complicated interactions of code blocks)
- Models can be built compositionally (we can naturally 'glue' models together at a high level, without worrying about implementation details and edge-cases)
- Models can be _migrated_ at a high level, whether interpersonally (collaboration with others who have a different ontology/vocabulary) or intrapersonally (one updates one's own model of the of world and wishes to reuse one's old model under new assumptions, without having to manually refactor code and dig into implementation details).

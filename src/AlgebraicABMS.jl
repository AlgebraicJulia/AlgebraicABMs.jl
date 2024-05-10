""" Some description of ths package
"""
module AlgebraicABMS

export hello

using Catlab

""" hello(name::String)

Returns the string "Hello, <name>!" where `<name>` is replaced with the provided parameter
"""
hello(name::String) = string("Hello, ", name, "!")

end

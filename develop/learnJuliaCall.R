library(JuliaCall)
julia <- julia_setup()

julia_command("a = sqrt(2)"); julia_eval("a")
julia_eval("sqrt(2)")
julia_call("sqrt", 2)
julia_eval("sqrt")(2)
## You can use `julia_exists` as `exists` in R to test
## whether a function or name exists in Julia or not
julia_exists("sqrt")
## You can use `julia$help` to get help for Julia functions
julia_help("sqrt")
## You can install and use Julia packages through JuliaCall
#julia_install_package("Optim")
julia_install_package_if_needed("Optim")
julia_installed_package("Optim") # version
julia_library("Optim")


julia_eval("import(Pkg)")
julia_eval('Pkg.activate("d:/julia/sesam.jl")')
julia_library("sesam")
julia_library("LabelledArrays")
#julia_eval('using sesam')
julia_eval("sesam.greet()")
julia_call("sesam.greet")

julia_eval("const se = s3a")
x = julia_call("se.getStateVector")
names(x) = julia_eval("symbols(se.getStateVector())")

xR <- c(B = 5, c = 4)
julia_assign(x, xR)

julia_eval("StateVectorType = @SLVector Float64 (:B, :R, :RN, :L, :LN, :I, :α)")

julia_createStateVector = function(x, module, name = toString(substitute(x))){
  #name
  cmd = str_interp("${name} = ${module}.StateVectorType(${name})")
  julia_assign(name,x)
  julia_eval(cmd)
  cmd
}
(julia_createStateVector(x, "se", "xs"))

f1 <- JuliaCall::julia_eval("
function f(x, names )
  LVector()
end")

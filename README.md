![STORMM logo](images/StormmLogoSmall.jpg?raw=true)

This code base will provide accessible, performant, interoperable libraries in a new family of MD
programs.  Basic operations such as coordinate and topology intake, user input parsing, and 
energy evaluations are managed through a common set of C++ classes and CUDA or HIP kernels.  In
many instances, HPC applications will be constructed using the same data arrays and layouts as
their CPU counterparts, although the most efficient GPU code may require sacrifices in the
convenience or performance of CPU code.  These sacrifices are likely to be minor, however, as the
best CPU code takes advantage of similar vectorization as GPU code and, even if vectorization
cannot happen in the same manner, the contents and format of basic data arrays do not restrict the
actual algorithm in any severe way.

---------------------------------------------------------------------------------------------------
  Installation Instructions
---------------------------------------------------------------------------------------------------
Here are some basic instructions:

- Download the repository and unpack the source in ```/your/STORMM/source/dir/```
- Set environment variables:
```
    STORMM_SOURCE  set to /your/STORMM/source/dir/
    STORMM_HOME    set to /your/STORMM/source/dir/
    STORMM_BUILD   set to /your/STORMM/build/dir/
    STORMM_VERBOSE set to COMPACT for concise testing output or FULL for more detailed results if
                   test programs are run manually (this can help to trace errors)
```
- Have CMake version 3.18 or greater, for CUDA compatibility.
- Make a new directory ```/your/STORMM/build/dir/``` and go to it, then type
  ```cmake ${STORMM_SOURCE} -DCMAKE_BUILD_TYPE=RELEASE [ additional cmake variable definitions ]```
  **OR**
- Execute ```cmake -b ${STORMM_BUILD} -s ${STORMM_SOURCE} -DCMAKE_BUILD_TYPE=RELEASE [ additional scmake variable definitions ]```
- Go to the ```${STORMM_BUILD}``` directory, if you are not already there.
- Type ```make -j``` to build with all possible resources.
- Type ```make test``` if the compilation is successful to run the test suite

There are several additional definitions that one can apply to the build.  Each is provided with
the cmake command prefaced with ```-D``` ("define").  They include:

- ```-DCMAKE_BUILD_TYPE``` with values ```RELEASE``` or ```DEBUG```.  The ```DEBUG``` version will
  eliminate optimizations present in the ```RELEASE``` version and engage debugging flags (that is
  ```-g``` to a C++ or Fortran compiler).  We use ```valgrind``` a lot when tracing errors in
  STORMM development, and the debugging flags are essential to getting line numbers in functions
  that exhibit memory errors or other exceptions.
- ```-DSTORMM_ENABLE_CUDA``` with values ```ON``` or ```OFF``` (default ```OFF```).  This will
  toggle the STORMM installation between CUDA-enabled and CPU-only versions.  The CPU-only version
  implements most features of the CUDA-enabled version, but some of the algorithms and precision
  models are not as expansive.  The intent is that, in the long run, the CPU version will be able
  to replicate all of what the GPU does using equivalent methods, albeit more slowly.
- ```-DSTORMM_ENABLE_RDKIT``` with values ```ON``` or ```OFF``` (default ```OFF```).  This will
  enable RDKit to communicate with STORMM, but requires a valid RDKit installation and does not yet
  unlock any real functionality.  Some features of STORMM are a stepping stone to information
  processing with RDKit, and some features of RDKit will enhance the capabilities on the STORMM
  roadmap.
- ```-DCUSTOM_GPU_ARCH``` with many values separated by semicolons (```;```) with no spaces.  Apply
  whatever architectures will be necessary for the GPUs available.  ```52``` will serve Maxwell
  GPUs such as the GTX 980, ```61``` will serve consumer-grade Pascal GPUs such as the GTX 1080 Ti,
  ```75``` will serve Turing GPUs such as T4 and the common RTX 2080 Ti, ```86``` will serve
  consumer-grade GPUs in the Ampere line (RTX 30xx (Ti), A40), ```89``` will serve consumer-grade
  GPUs in the Lovelace line (RTX 40xx (Ti)), and settings such as ```60```, ```70```, or ```80```
  will serve X100 series data-center grade GPUs suchas GP100, V100 ("Volta"), and A100.  The Hopper
  line (H100, input value ```90```) is not well tested and some features were encountering errors.
  The default behavior, triggered in the absence of ```-DCUSTOM_GPU_ARCH``` specification, is to
  set the native ```CMAKE_CUDA_ARCHITECTURES``` variable to ```52;60;61;70;75;80;86;89```.
  Trimming the list, particularly with regard to older GPU hardware, will reduce build times.  The
  new variable ```CUSTOM_GPU_ARCH``` is provided because the native ```CMAKE_CUDA_ARCHITECTURES```
  comes with a preset of ```52``` which will omit a lot of optimizations on more recent cards.
- ```-DCMAKE_SHOW_PTXAS``` with values ```ON``` or ```OFF``` and default ```OFF```.  This will
  cause the compiler to output ```--ptxas``` output to the terminal, which can be helpful in
  gauging kernel performance (in particular, regarding register usage).  Setting this to ```ON```
  will output 3-4 lines for every version of every kernel compiled, which can be a lot to parse as
  some of the lines are 240+ characters in length.
- ```-DSTORMM_ENABLE_TEST_COVERAGE``` which is undefined by default but can be defined with a valid
  installation of ```gcov``` to check the coverage of unit testing code during development.

---------------------------------------------------------------------------------------------------
  Code Standards
---------------------------------------------------------------------------------------------------
- This code base will follow the C++17 standard, in particular for enum classes, member
  initializers, brace initialization for POD aggregates, deleted and defaulted functions, nullptr,
  constructor delegation, and key Standard Library algorithms.  C++17 is comprehensible and
  compatible to NVIDIA's compilers.  It subsumes the major C++11 increment in C++, will be
  insulated from deprecation for years to come, and remains a safe choice for hybrid programming.
- Developers should strive for a code feel much like C, with basic C++ features like templates,
  private or protected struct or class member variables with public getter functions, and often
  passing by reference rather than by pointer (in C++, type & rather than * in C).  The code base
  should be familiar to a developer who understand functions, loops, structs, and pointers.  Such
  a developer will need to learn (struct or class) member functions, enumerator class types,
  initialization of objects, templates, and function overloading.  These concepts can be learned
  independently and may reinforce one another.  The transition from intermediate C or python coding
  to writing C++ with STORMM is intended to be a steady learning curve.
- Iterators, ranged for loops, and unique pointers are among the features of modern C++ deemed to
  be further from the original C language than most programmers should have to go, and therefore
  seldom appear in the STORMM libraries themselves.  Programmers writing backend or small programs
  that link to the STORMM libraries may use these features at their discretion.
- All dynamic memory intended to operate in a hybrid CPU / GPU environment should be allocated with
  Hybrid objects from the STORMM library.  This class functions with many of the features of the
  STL vector container and interfaces with it for easy conversion to any function requiring that
  type.

---------------------------------------------------------------------------------------------------
  Coding Conventions
---------------------------------------------------------------------------------------------------
- We're not doing this with a formatter!  The goal is to get as much code on one screen as possible
  in a clear and legible fashion.  Code formatters tend to insert excessive line breaks where they
  are not necessary or make it harder to skip lines where it would be nice to separate groups of
  statements.  They also work against the concept of column uniformity in assignment operators,
  which can be helpful within a scope or a tightly connected group of statements.
- Max line width: 99 characters (100th must be a line break).  This may be violated if a kernel
  launch cannot fit within the space, but in most cases names should be adjusted to avoid wrapping
  past the 100 character limit.
- Indent new scopes by two spaces.
- When continuing arithmetic statements on new lines, end each line with an arithmetic operator,
  which stands out as an indication that something must follow.  Indent the next line as far as the
  assignment operator, i.e.:

```
    t  = x + y +
         z + w;
    t += (k * p) +
         (x / y) + 4;
```

- When writing arithmetic statements, add parentheses and group operations according to the order
  of operations.  See the example above.  Use spaces between each arithmetic operator and its
  operands, but do not add extra spaces between variables and left or right parentheses.
- When making arithmetic statements within index assignments, or with very short NON-STRUCT
  variable names, it is permissible to show the order of operations by cutting the spaces away from
  basic arithmetic operators.
- When continuing lines within parentheses, brackets, or braces, indent to the interior of the most
  recent open delimiter.
- Align assignment operators for a tightly connected series of statements, as shown in the code
  above.
- Use of auto-typing (type inference) is discouraged, to make the code as inviting and clear as
  possible for new developers.
- When writing conditional statements that branch the code path when a logical expression
  logic_expr is false, use (!logic_expr) when performance is critical.  Use (logic_expr == false)
  in most other cases: this implies an extra comparison, whereas ! is merely arithmetic, but the
  (logic_expr == false) syntax is more legible when the meaning of the code is most important.
- When passing by reference or declaring a reference to an object, attach the & to the left side of
  the object name, rather than to the right side of the object name, i.e. int foo(MyClassType &mc)
  not int foo(MyClassType& mc) or int foo(MyClassType & mc).  The middle choice can let the
  reference symbol get lost, while the final choice can be confused with the bitwise AND operator.
- When returning a reference from a function, attach the reference to the right side of the
  returned type, not to the left side of the function name.
- When printing particular formats of numbers, i.e. 0.0F or 7LL, use capital letters in particular
  to prevent lowercase "l" from being confused with one (1).

---------------------------------------------------------------------------------------------------
  Abstraction and Dependencies
---------------------------------------------------------------------------------------------------
- Class inheritance, while powerful, has not yet been implemented for STORMM objects.  The most
  useful and appropriate places would likely be among the various coordinate objects, but even
  there the practice would be of limited benefit as each of the differet coordinate objects are
  somewhat orthogonal.  In general, STORMM development does not prohibit inheritance but does not
  require it, either.
- For including libraries, follow this order:
  * C++ Standard Libraries (```#include <vector>```), in alphabetical order
  * ```#include "copyright.h"``` (the copyright statement is found in the base source directory)
  * STORMM libraries (```#include "Math/vector_ops.h"```), in alphabetical order
- C++ implementation (.cpp) files should generally include their own corresponding header files.
  Any classes, functions, or other objects present in a header file should be traceable by some
  ```#include``` statement in that header file.  Any classes, funcions, or other objects present
  in an implementation file should be traceable by additional ```#include``` statements in the
  implementation file itself, or in the relevant library header file.
- Place documentation for function usage in header (.h) files and implementation-relevant
  documentation in implementation (.cpp or .tpp) files.  Place the implementations of template
  functions in template implementation (.tpp) files, not in the header file where the templated
  function or class template is first declared.  The .tpp file should be included at the end of the
  corresponding header file, outside of any namespace scopes but within the header file's customary
  ```#include``` guards.  Template files should re-declare the relevant namespace scopes to wrap
  their implementations, but not contain additional ```#include``` guards.

---------------------------------------------------------------------------------------------------
  Naming Conventions
---------------------------------------------------------------------------------------------------
- Camel case (thisIsCamelCase) for function names
- Lowercase lettering with underscore separators (these_are_underscores) for variable names
- Uppercase lettering with underscore separators (THIS_IS_A_MACRO) for macro functions and
  macro definitions
- Capitalization of first letters (ThisCapitalizesAllFirstletters) for class names
- Use all uppercase letters for enumerations, with underscore separators (THIS_IS_AN_ENUMERATION).
  Avoid #define'd constants, preferring constexpr, but when necessary maintain the enumeration
  case convention for #define'd constants.
- Keep variable names descriptive, but not excessively verbose.  Names of major classes and
  critical variables should contain complete words, avoiding abbreviations except for obvious
  diminiutions like "Information" -> "Info" or "Matrix" -> "Mat".
- The length of the variable name should follow from the extent of the scope.  Short variable names
  are permissible in very local contexts, longer ones necessary for expansive scopes.
- File names should be in all lowercase with underscore separators, similar to variable names.
- Library files get .cpp extensions, header files .h extensions.  Do not use separate header files
  for data structures and functions within a library--if circular referencing occurs, make a new
  library and separate out the problematic data structures.

---------------------------------------------------------------------------------------------------
  Function Declarations
---------------------------------------------------------------------------------------------------
- Pass by class objects by const reference to C++ functions, unless they will undergo changes in
  the function--then, pass them by pointer to indicate that they may change as the function is
  executed (see the discussion on pointer usage below).  Pass abstracts (structs of const pointers
  and scalars) by value to CUDA or HIP kernels.  Depending on the situation, passing an abstract by
  value to a C++ function (to avoid an extra layer of de-referencing) is permissible.  
- Place function descriptions in the .cpp library file, ahead of the function declaration, not in
  the header.
- Declare const-ness of elementary type variables in the library file but not in the header file,
  i.e. const double x in the .cpp but double x in the .h.  Exceptions to this rule include template
  header files, where there may be no corresponding code in a .cpp file.
- Place doxygen \brief comments for classes and structs in the .h file above their declarations.
  Place \brief comments for functions or extern / global class objects above their code in the .cpp
  library file.
- Declare overloaded functions using comments such as

```
    /// \brief General description of function in all guises
    ///
    /// Overloaded:
    ///   - First version
    ///   - Second version
    ///
    /// \param (Describe all parameters to any of the function's overloaded versions)
```

- Overloaded functions should cascade from instance to instance so long as it does not create
  ineffciency in production calculations.  Conserve code and avoid replication, but not at the
  expense of performance.  Use the "buck stops here" overloaded function in the most performance-
  critical sections of algorithms for versatility as well as performance.

---------------------------------------------------------------------------------------------------
  Struct Declarations
---------------------------------------------------------------------------------------------------
- Place descriptions for structs in the .h header file where the struct is declared.
- Most structs will have private member variables that are accessible only by getter and setter
  functions.  In such cases the developer will often define their own constructor.  For
  convenience, such as with a small struct which serves only to interlace related data in an array
  or pass by value to another function, structs may be "unguarded" meaning that they have all
  public member variables.  Such structs should also have no member functions--they are, for all
  practical purposes, structs in the traditional C sense, but this is of course a convention, not
  a compiler-enforced rule.
- Template structs take template type name T, or T1, T2, T3... as necessary.  When multiple
  templated types are present, tuples of the base type should be named with the size of the tuple,
  e.g. ```<typename T, typename T4>``` when T could be ```int``` or ```float```.  More elaborate
  names are encouraged when the types serve different purposes and are not expected to comprise
  similar number formats, e.g. ```<typename Tcoord, typename Tcalc>``` for distinguishing the data
  types used to express coordinates and the data type used to perform most arithmetic.
- Use initializer list mechanics where possible, but for structs that require detailed reading of
  a text file or other complex manipulation of data, write an argumentless constructor that builds
  a blank struct, then use that to fulfill the initializer lists of any other constructors which
  take actual arguments (i.e. a file name or other data struct).  Use the function bodies of the
  subsequent constructors to resize arrays and other member variables as necessary and fill out
  the contents of the struct.
- Use cascading initalization, including the above strategy, as much as possible.  Assume that
  struct initialization is not a performance issue, and write code such that it is an overhead
  cost but not a factor in calculation production.

---------------------------------------------------------------------------------------------------
  Template Conventions
---------------------------------------------------------------------------------------------------
One of the most significant features that distinguishes C++ from the original C language is the
existence of templates, objects or functions which are (within limits) agnostic to the data type
that they will operate upon.  Because they step into the realm of modern C++ and a new world of
abstraction and inference, templates make an important "boundary" case of what can be sanctioned
within the (arbitrary, self-imposed) STORMM coding conventions.  The critical aspect of successful
template design concerns type inference, a depthy and general aspect of C++ which made laudable
improvements in C++11.  When a template is _used_ in the code, the compiler essentially creates an
overloaded variant of the template function as stipulated either by an explicit type declaration,
i.e. ```double t = foo<double>(... argument list ...)``` or by infering the exact overloaded type
from a function argument of the template type, i.e. a primitive such as
  
```
    template <typename T> T foo(const std::vector<T> &list) {
      ... function content ...
    }
```

In the above case, the return type of ```foo()``` cannot be dictated by the type of variable
that _receives_ the output.  However, it can be inferred from the argument fed to ```foo()```:

```
    std::vector<long long int> phone_numbers(1000, 0);
    ... load phone_numbers with data... 
    float local_var = foo(phone_numbers);
```

The template leads to an overloaded instantiation of the function that takes a vector of ints,
and therefore returns an int, regardless of the type of ```local_var```.  This can be very
powerful, but it can also lead to some mistakes that will break code later on.  Let's say that
```foo()``` was designed to accept a vector of int4 objects:

```
    struct int4 {
      int x;
      int y;
      int z;
      int w;
    };
```

If there were any code inside of ```foo()``` that operated on components of the input data such
as x, y, or w, this would negate the use of ```foo()``` on any data type that does not have those
components.  Similarly, if ```foo()``` contained an operation like ```+=``` that was well defined
for standard ints but not for int4s, its use would be diminished.  All of these errors would be
silent until the templated function was given a type to work on, at which point it would either
compile or not.  With this understanding:
- Write templates such that they can be used with inputs of scalar data, or vectors of scalar
  types, with type inference by the compiler rather than explicit declaration.
- For templated functions that take no input arguments, explicit declaration of the intended
  type is required as a function modifier, i.e. ```foo<int>()```.  Write these functions to be
  able to work with many types, without operations that would restrict the input type, or with
  names that suggest the types that the function is intended to work with.
- Understand that templates cascade when one function call another, and that this behavior
  cannot be protected by an ```if (conditions) { }``` statement in the code.  Any code that
  appears in a templated function must be able to work with any type that might be presented to
  the template.
- For templates that are only intended to work with a finite number of types, use functions with
  traps to throw runtime errors when developers attempt to instantiate the template with the
  wrong type.  This feature will be superceded by code that throws compile-time errors in the
  future.

---------------------------------------------------------------------------------------------------
  References and Pointers
---------------------------------------------------------------------------------------------------
In C, it is common to pass arrays and other large objects by pointer.  This saves the program from
having to create a copy of the argument in question solely for the purposes of that function, at
the expense of de-referencing the pointer at every access of the object or one of its attributes.
This also requires frequent use of the -> operator, i.e. ```x->y``` for member variable y in
struct x, if x has been passed by pointer (```foo(*x)```).  In C++, it is more common to pass by
reference (```&x```).  There are some subtle differences between references and pointers, but from
the standpoint of performance references and pointers are equivalent (the program still needs to
de-reference the reference).  Passing by reference also removes the need to write -> for accessing
object attributes and also &x in the function call (it can be confusing, as to pass by reference
there will be a &x in the function declaration).  One can write ```x.y``` inside a function passed
argument ```&x``` (although the de-referencing will still occur), and call foo(x).  It
is not best _just_ to pass arguments by reference, however: whenever passing by reference, STORMM
passes by const reference, restricting the function from changing the data referenced by the
variable (the only thing that a function can produce is its return value).  While some developers
strive to have every argument be const (including arguments passed by value), there are situations
where it is far preferable to have a function do two things, usually modifying some array while
accumulating a result.  In these cases, the best practice is to return the result and pass the
variable which will be modified through the process by pointer.  That way, in the function call,
developers will need to write ```&``` when passing the argument, sending a clear signal that it,
too, will be modified over the course of the call.  For example:

```
(... Declare BuildingDesign and ConstructionProject objects ...)

//-------------------------------------------------------------------------------------------------
int squareFootage(const BuildingDesign &blueprint, std::vector<ConstructionProject> *buildings,
                  const std::vector<int> &addresses) {
  int sq_ft = 0;
  for (size_t i = 0; i < addresses->size(); i++) {
    buildings->data()[addresses[i]].floor_count = blueprint.floors;
    buildings->data()[addresses[i]].entry_dimensions = blueprint.door_width;
    buildings->data()[addresses[i]].length = blueprint.length;
    buildings->data()[addresses[i]].width = blueprint.width;
    sq_ft += (blueprint.floors * blueprint.length * blueprint.width) - (blueprint.door_width * 2);
  }
  
  return sq_ft;
}

//-------------------------------------------------------------------------------------------------
int main() {
  BuildingDesign home_blueprint1;
  BuildingDesign home_blueprint2;
  std::vector<ConstructionProject> my_town;
  std::vector<int> home_addresses1;
  std::vector<int> home_addresses2;
  (... Fill in the details of home_blueprint1, home_blueprint2, and home_addresses(1,2) ...)
  
  // The function call below will accumulate the total residential space while assigning specific
  // addresses to be a particular type of home.
  int total_residential_space  = squareFootage(home_blueprint1, &my_town, home_addresses1);
  total_residential_space     += squareFootage(home_blueprint2, &my_town, home_addresses2);
  return total_residential_space;
}
```

The program above passes blueprint objects and the arrays of home addresses by const reference,
while passing the array of actual buildings in the town by pointer.  In the main function, it is
obvious that ```my_town``` is being modified as the details of particular buildings are filled out
by calls to a function that also accumulates the square footage of buildings made with a certain
blueprint.  It would also be feasible to write the ```squareFootage()``` function with a non-const
reference for its second argument, but then the calls in main would give no indication that
```my_town``` was being modified by each call: const and non-const references look the same when
calling a function, but passing by pointer guards against stealth data modifications.

Do not be afraid of pointers!  C++ is a superset of C, and therefore can make use of them.
It is not desirable to use pointers to dynamically allocate memory, however: that is for
containers and smart pointers, the most common being ```std::vector```.  For array access, pointers
are often (but not always) less desirable than containers due to the fact that a pointer to a block
of memory is just that: it carries no information about what the bounds are, and while it is
possible to reach a segmentation fault by accessing a non-existent element of a std::vector it is
much easier to go into the code and install a bounds check if that happens.  In optimized code,
accessing an element of a std::vector by the ```[]``` operator is as fast as accessing an element 
of an array through a pointer--the compiler knows what the program is meant to do, so pointers
seldom have a performance advantage unless they help to skip over an extra layer of function calls
and de-referencing (there are cases of that, and STORMM libraries frequently provide small structs
full of const pointers to expedite access to data, especially in GPU-based arrays).  However, for
the specific purpose of letting a function modify one of its input arguments, pointers are the
preferred route.

---------------------------------------------------------------------------------------------------
  Class Abstracts and the C++ >> C >> HPC Transition Model
---------------------------------------------------------------------------------------------------
STORMM is a collection of C++ classes and CUDA (and perhaps in the future HIP or OpenCL) Kernels
which communicate by stripping the data down to a lowest common denominator of C.  Each class
built in STORMM should emit one or more "abstracts" which will collect critical constants, sizing
information for arrays, and raw pointers to the underlying information contained in the class.
Each abstract is a ```struct``` and as such should contain only a basic constructor taking values
for each member variable as a formal argument, copy and move constructors, and perhaps a manually
defined destructor.  This represents a C++ to C transition: it is possible to overrun the limits
of various arrays when uing an abstract if the stated bounds are not respected, as there are no
accessor functions with their own internal bounds checks as would be expected in a C++
```class```.  Abstracts should contain suffixes like ```Reader```, ```Writer```, ```Kit```, or
```Plan``` depending on what seems most appropriate given the nature of the object and the
mutability of data mapped by the abstract.  In most cases, even ```Writer``` abstracts will contain
some ```const``` members such as array sizing constants, and thus will not be able to have copy and
move assignment operators.  In most cases the abstracts will be obtained from class member
functions named ```data()``` but more descriptive names may be used to produce specialized
abstracts tailored to a subset of the underlying class object's contents.

---------------------------------------------------------------------------------------------------
  Templated CUDA Kernels and Templated Kernel Arguments 
---------------------------------------------------------------------------------------------------
A thought experiment can establish that it is not feasible to call a templated CUDA kernel from a
C++ file (.cpp) unless the HPC compiler is also compiling the .cpp implementations, i.e. is the
only compiler in town.  This is not the way STORMM is intended to function, given the standalone
CPU functionality that most of the code base aims to provide for prototyping and basic testing.
This does not, however, mean that a templated class instantiated by a C++ compiler cannot be
delivered to a CUDA (or HIP, or OpenCL) implementation without writing separate branches for every
intended data type of the template, an approach which would be prohibitive in the case of multiple
independent templated data types.  The preferred approach in STORMM is to create abstracts of the
templated class with the templated data pointers cast to ```void*```.   This cannot be done with
scalar constants included in the abstract, but it is seldom a problem.  Abstracts of this sort are
typically obtained by member functions named ```templateFreeData()```.  The void-casted abstracts
are then converted back to their original types by templated free functions named
```restoreType```, overloaded by the different abstracts they take for formal arguments.  Because
both C++ and HPC compilers will build objects out of the various header files, each will develop
their own template instantiations for a given object while the ```void*```-casted form is useful
as the intermediary.

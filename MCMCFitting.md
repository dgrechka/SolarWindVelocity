Writing F# code
---------------
In standard Markdown, you can include code snippets by 
writing a block indented by four spaces and the code 
snippet will be turned into a `<pre>` element. If you do 
the same using Literate F# tool, the code is turned into
a nicely formatted F# snippet:

    /// The Hello World of functional languages!
    let rec factorial x = 
      if x = 0 then 1 
      else x * (factorial (x - 1))

    let f10 = factorial 10

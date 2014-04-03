options(rstudio.markdownToHTML = 
          function(inputFile, outputFile) {      
            require(markdown)
            markdownToHTML(inputFile, outputFile, stylesheet='github.css')   
          }
)

# source: https://gist.github.com/andyferra/2554919
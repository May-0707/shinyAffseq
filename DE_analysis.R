### when platform changes, update clinical table
observe({
  #shinycat("observe platform to update clinical table...\n")
  if (is.null(input$GEOaccID)) {
    ex = NULL
    clinMatrix = NULL
    return(NULL)
  }

})

######################################
##  different groups
######################################

###################################################
# get possible values of the selected column names
###################################################

groupsForSelectedColumn <- reactive({
  #shinycat("In groupsForSelectedColumn reactive...\n")
  vars = clinMatrix #values.edit$table
  if (is.null(vars) | is.null(input$Groups)) {
    return(NULL)
  }
  
  vars <- vars[, as.character(input$Groups)] 
  vars = factor(vars)
  return(as.list(levels(vars)))
})


#observeEvent(input$groupGO, ({
#  #closeAlert(session, "merge-alert")
#  content = "Merge two or more groups together by selecting the groups from the drop-down boxes on the left, and specifying a name for the new group in the corresponding text boxes on the right. Click on Save, and a new column will be added to the clinical data table."
#  add = paste0("<p><p> Selected column: <strong>", input$selectedColumn, "</strong></p>")
#  
#  content = paste0(content, add)
#  
#  #createAlert(session, "mergeGroupsAlert", alertId = "merge-alert", title = "Current Status: Merging", style = "shinygeo-primary", content = content)
#})
#)

#output$selectedGroups <- renderUI({
#  selectInput('newGroups','Select or Merge datas to Two Groups for Comparison', 
#              choices = groupsForSelectedColumn(), multiple=TRUE,
#              selected = defaultGroupsForSelectedColumn(),
#              width='100%',
#              selectize = TRUE
#              
#  )
#})
####################################################################
## renders drop-down menus (server-side) for clinical group 
## selection for merging groups in MergeGroupsModal
####################################################################

#output$groupA <- renderUI({
#  selectInput('selectGroup1', 'Group 1(e.g. Control)', 
#              choices = groupsForSelectedColumn(), multiple=TRUE,
#              selected = NULL,
#              width='80%',
#              selectize = TRUE
#  )
#})
#
#output$groupB <- renderUI({
#  selectInput('selectGroup2', 'Group 2', 
#              choices = groupsForSelectedColumn(), multiple=TRUE,
#              selected = NULL,
#              width='80%',
#              selectize = TRUE
#  )
#})


####################################################################
## respond when Save button is clicked on MergeGroups modal
####################################################################
#observeEvent(input$GetDEGs, ({
  #shinycat("Merging groups...\n")
  #content = ""
  #col = "selfDefGroup"
  #
  #g1 = input$selectGroup1
  #g2 = input$selectGroup2 
  #
  #g1 = g1[g1!=""]
  #g2 = g2[g2!=""]
  #
  #g.all = c(g1,g2)
  
#  if (length(g.all) > length(unique(g.all))) {
#    content = paste0(content, 
#                     "<p> Error: A value cannot appear in multiple groups <p>")
#  }
#  
#  if (col%in%colnames(clinMatrix)) {
#    content = paste0(content, 
#                     "<p> Error: Column Name Exists. Please select a new column name <p>")
#  }
#  
#  if (content!= "") {
#    createAlert(session, "mergeGroupsAlert", alertId = "merge-alert-error", title = "Save Error", style = "shinygeo-danger", content = content, append = TRUE)
#    return(NULL)
#  }
  
  #clin[which(clin[,groupcol] %in% c(g1)),selfDefGroup] = "Group1"
  #clin[which(clin[,groupcol] %in% c(g2)),selfDefGroup] = "Group2"
  
  #X = as.character(values.edit$table[[input$selectedColumn]])
  #Y = rep("", length(X))
  #
  #add1 = "## merge groups from selected column ##\n"
  #add1 = paste0(add1, "tmp = as.character(data.p[[\"", input$selectedColumn, "\"]])\n")
  #add1 = paste0(add1, "Y = rep(\"\", length(tmp))\n") 
  #
  #if (length(g1) > 0 & input$group1Label != "") {
  #  Y[X %in% g1] = input$group1Label 
  #  add1 = paste0(add1, "Y[tmp %in% ", vector.it(g1), "] = \"", input$group1Label, "\"\n")  
  #}
  #if (length(g2) > 0 & input$group2Label != "") { 
  #  Y[X %in% g2] = input$group2Label 
  #  add1 = paste0(add1, "Y[tmp %in% ", vector.it(g2), "] = \"", input$group2Label, "\"\n")  
  #}
  #
  #if(length(unique(Y)) <= 1) {
  #  createAlert(session, "mergeGroupsAlert", alertId = "merge-alert-error", title = "Save Error", style = "shinygeo-danger", content = "<p> Error: this merge would create a column where all values are the same. This operation is currently not supported" , append = TRUE)
  #  return(NULL)
  #} 
  #
  #data = values.edit$table
  #
  #data[[col]] = Y
  #
  #
  #add1 = paste0(add1, "data.p[[\"", col, "\"]] = Y\n")
  #isolate(add.code(add1))
  
  #isolate(values.edit$table <- data)
  #toggleModal(session, "MergeGroupsModal", "close")
  #updateSelectInput(session, "selectedColumn", choices = ColumnNames(),
  #                  selected = col)
#})
#)

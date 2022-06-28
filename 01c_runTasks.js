// 1. Execute o script para gerar os tasks do jeito exato que você quer que sejam rodados (escala, etc)
// 2. Quanfo todas as tarefas estejam com o ícone ´Run´ em azul na aba de task, aperte F12 para entrar no modo de desenvolvedor do Chrome
// 3. Na opção ' Console' digite:
runTasks1= function(){
    $$('.run-button' ,$$('ee-task-pane')[0].shadowRoot).forEach(function(e) {
         e.click();
    })
}
runTasks2= function(){
    $$('ee-table-config-dialog, ee-image-config-dialog').forEach(function(e) {
    var eeDialog = $$('ee-dialog', e.shadowRoot)[0]
var paperDialog = $$('paper-dialog', eeDialog.shadowRoot)[0]
    $$('.ok-button', paperDialog)[0].click()
    })
}
// 4. depois digite:
runTasks1();
// 5. e quando apareçam todas as janelas de confirmação digite:
runTasks2();

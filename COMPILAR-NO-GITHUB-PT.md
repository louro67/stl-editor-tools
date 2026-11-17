# Compilar o STL Editor Tools 0.3.1 no GitHub (sem instalar Rust no PC)

1. Entre na sua conta em https://github.com.
2. Crie um repositório novo, por exemplo `stl-editor-tools` (pode ser privado).
3. Extraia este ZIP no computador.
4. No repositório novo, escolha **Add file > Upload files** e arraste o conteúdo da pasta `mesh-solidify-0.2.0` (incluindo a pasta oculta `.github`). Não carregue apenas a pasta ZIP.
5. Clique em **Commit changes**.
6. Abra o separador **Actions** do repositório.
7. Escolha **Build STL Editor Tools** na lista à esquerda e clique em **Run workflow**.
8. Aguarde que a execução termine. Verde significa que compilou; vermelho significa que ocorreu um erro.
9. Abra a execução concluída. Na parte inferior, em **Artifacts**, descarregue `STL-Editor-Tools-0.3.1-WASM`.

O ZIP do Artifact é o pacote que o script de compilação produziu e validou com `wasm-tools`. Se a execução falhar, abra a execução, copie o primeiro erro vermelho (ou envie uma captura de ecrã) para podermos corrigir o código.

**Nota:** este processo compila e valida o componente WASM; não garante por si só que todos os comandos sejam aceites pela instalação específica do SolidCraft 0.43. É necessário testar no SolidCraft.

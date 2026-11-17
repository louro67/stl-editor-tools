# STL Editor Tools 0.3.1 — alterações ao código-fonte

Esta versão parte do código-fonte original `mesh-solidify 0.2.0` e acrescenta entradas que encaminham operações para comandos nativos do SolidCraft.

## Novas operações no menu Plugins

- **Convert Mesh to Solid** — chama `mesh.convert` (malha fechada; tolerância 0,05 mm).
- **Scale Body** — chama `solid.scale`; escala uniforme ou independente nos eixos X/Y/Z.
- **Split Body by Plane** — chama `solid.split_body`; planos XY, XZ ou YZ e deslocamento em mm.
- **Delete Face** — chama `face.delete`; indicar um ponto XYZ situado na face; opção de reparar/fechar a abertura.
- **Press Pull Face** — chama `solid.press_pull`; indicar um ponto XYZ na face e uma distância em mm.
- **Move Body** — chama `solid.move`; deslocamento X/Y/Z em mm.

Mantém **Solidify Mesh** e **Inspect Mesh** originais.

## Limitações importantes

Estas operações usam o motor geométrico nativo do SolidCraft. Não são um editor de vértices/arestas em modo livre. `Delete Face` e `Press Pull Face` identificam uma face por um ponto XYZ, pelo que é necessário saber um ponto sobre a face. `mesh.convert` requer uma malha fechada. A API do plugin não fornece acesso livre para substituir diretamente a malha arbitrariamente.

## Compilar o plugin instalável

O ZIP contém o código-fonte alterado e o manifesto, mas não inclui um `.wasm` recompilado. Para gerar o plugin instalável é necessário:

1. Rust 1.85 ou superior e o target `wasm32-unknown-unknown`.
2. `wasm-tools` instalado e disponível no PATH.
3. Executar `scripts/build-plugin.sh` na raiz do projeto.

O resultado esperado é `dist/mesh-solidify-plugin.zip`. Não instalar este ZIP de código-fonte no SolidCraft: primeiro é necessário compilar.

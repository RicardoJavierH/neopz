// -*- c++ -*-

//$Id: meshes.h,v 1.2 2005-12-19 12:19:36 tiago Exp $

class TPZCompMesh;
#include "pzvec.h"

/** Barra tracionada - Validado
 */
TPZCompMesh * BarraTracionada(int h, int p); 

/** Barra tracionada girada de 10 graus em torno do eixo x.
 * int cubo = 1 faz um cubo - Validado
 * int cubo = 2 faz um prisma - Validado
 */
TPZCompMesh * BarraTracionadaGirada(int h, int p);
 
/** Viga engastada com momento concentrado na extremidade livre - Validado
 */
TPZCompMesh * VigaEngastada(int h, int p);
void MomentoExtremidade(TPZVec<REAL> &pto, TPZVec<REAL> &force);

/** Viga engastada com forca de corpo (volumetrica) - Validado
 */
TPZCompMesh * VigaEngastadaForcaVolume(int h, int p);


 
/** Viga bi-engastada
 */
// TPZCompMesh * VigaBiEngastada(int h, int p);
 
// /** Vigas em balanćo
//  */
// TPZCompMesh * VigasBalanco(int h, int p);
//  
// /** Bloco confinado
//  */
// TPZCompMesh * BlocoConfinado(int h, int p);

void SetPOrder(int p);

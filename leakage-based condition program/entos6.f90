!----------------------------------------------------------------------------------------------------------------------------------------------------
!MHICO - Modelo HIdrodin�mico COnjunto.
!*** Versi�n Mhico3 con cambios (denotados por !***, las l�neas comentadas as� y sin explicaciones son l�neas a�adidas) respecto a Mhico
!Programa MEF v�lido para mallas de hasta 99999 elementos o 99999 nodos (25/02/2013). Autor: H�ctor Garc�a R�bade
!Lectura de ficheros .txt (conectividades, coordenadas, propiedades y CC) y generaci�n de ficheros .dat con formato por Tecplot	(Tecplot Data Loader)
!----------------------------------------------------------------------------------------------------------------------------------------------------

!Posibilidades del programa:
!---------------------------
!1. Se pueden resolver tres ecuaciones 2D diferentes: ec. NS2D (en funci�n del calado en vez de la presi�n, ec. te�rica sin conservaci�n de masa en 3D), 
! ec. aguas someras (integraci�n de ec. NS3D en altura, prop de manning) y ec. agua subterr�nea (integraci�n de ec. continuidad en altura, propiedades 
! de conductividad, porosidad y �ngulo de anisotrop�a) para acu�feros no confinados.
! En ec. NS2D se calculan las velocidades y el calado (igual a la altura de la l�mina), en ec. aguas someras se calculan velocidades y la altura de la 
! l�mina (el calado postproceso), en ec. agua subterr�nea se calcula la altura de la l�mina o nivel fre�tico (el espesor fre�tico y las velocidades 
! postproceso).
!2. La ec. aguas someras incorpora un t�rmino de fricci�n dependiente del n=n�mero de manning (sol. en canales con Re bajo y n=0 equivalente a sol. 
! con Re alto y cierto n).
!3. Uso de elementos de tres nodos en la ec subterr�nea y de elementos mixtos P2-P1 en la ec. NS2D y la ec. aguas someras.
!4. Tiene tres modelos de elementos finitos: superficial (ec. NS o aguas someras), subterr�neo (ec. subt) y conjunto del tipo superficial-subterr�neo (ec. 
! aguas someras y subt) desconectado (pero no colgado) para acu�feros libres.  !*** "desconectado (pero no colgado)" en vez de "conectado"
!5. Elementos interfaz - El modelo de aguas someras utiliza la condici�n seco-mojado, y el modelo conjunto utiliza una condici�n similar que utiliza CC 
! derivadas del intercambio de informaci�n entre modelos.
!6. Modelo conjunto con conservaci�n de masa en el contorno m�vil (uso de 1 s�lo contorno m�vil ya que lo permite, a�n en caso de trabajar con espesores 
! diferentes al haber el mismo n�mero de nodos).  !*** Se usan espesores diferentes en el contorno m�vil a los usados en Mhico
! C�lculo para cada modelo de alturas, espesores (calado o espesor fre�tico) y velocidades. En el modelo conjunto se representan conjuntamente 
! velocidades, alturas e individualmente el calado (valor nulo si hay espesor fre�tico) y el espesor fre�tico (valor nulo si hay calado). En el modelo 
! superficial se representan adem�s la vorticidad y las tensiones.
!7. Si se simula conjunto siempre se puede simular superficial (condici�n seco-mojado).
!8. Posibilidad de resolver un esquema impl�cito (para todas las ec; resoluci�n estacionaria conlleva a un problema no lineal, resoluci�n transitoria 
! conlleva un problema no lineal por incremento de tiempo) o uno semi-impl�cito de crank-nickolson (s�lo para modelo superficial, ec. NS2D o aguas 
! someras; resoluci�n siempre transitoria, que conlleva a un problema lineal por incremento de tiempo -soluci�n estacionaria si CC ctes-) 
!9. Posibilidad de usar el m�todo de newton (s�lo para ec. aguas someras). 
!10. C�lculo de las integrales de contorno de las tensiones viscosas (muchas veces despreciadas) para ec. NS2D o ec aguas someras
!11. Estabilizaci�n para ec. NS2D o aguas someras. Mejora la resoluci�n del sistema lineal y la convergencia del problema no lineal.
!12. Uso de Gradientes Biconjugados con un precondicionador diagonal (usando la norma de los coeficientes de la fila, PCGB) o con un precondicionador LU 
! aproximado (PCGBLU). Cambiando el c�digo (partes comentadas) tambi�n podr�a utilizarse PCGB sin precondicionamiento o con precondicionamiento con los 
! coeficientes de la diagonal.
!13. Se puede aplicar lluvia (ec. aguas someras y agua subterr�nea) y bombeos (ec. subt) por pantalla. 
!14. Los archivos .dat tienen el formato adecuado (en doble precisi�n) para que Tecplot pueda hacer una pel�cula con la soluci�n transitoria.
!15. Para el mismo n� de elementos el modelo superficial es m�s lento por haber mayor n� de ecuaciones y por uso de elementos P2P1. Incremento de la 
! velocidad de computaci�n con: 
! Utilizaci�n de vectores para el almacenamiento de la matriz dispersa del sistema (necesario en mallas de muchos elementos para tener espacio para 
! almacenar la matriz).
! Optimizaci�n en el almacenamiento de los coeficientes, mediante escritura de coeficientes en misma posici�n y su posterior reordenamiento de los vectores
! Uso de la �ltima soluci�n de Picard como aproximaci�n para resolver el siguiente lineal.
! Modelo superficial: posible uso para Re bajos uso del m�todo de newton (aguas someras) o uso de PCGBLU, posible uso para Re altos de estabilizaci�n.  
! Modelo subterr�neo: posible uso de PCGBLU
! Modelo conjunto: uso de 1 it sup por cada sol del modelo subterr�neo lo que es necesario porque permite acoplamiento y lo se�alado para cada modelo 
! (si se usa un m�todo de PCGB ser� para ambos modelos). 
!16. C�digo compatible con .f95. Ello permite su uso en compiladores de 64 bits con los que la ejecuci�n del programa ser� mucho m�s veloz.

!Motivos de error (datos de etiquetas para f�cil localizaci�n mediante b�squeda de etiquetas):
!---------------------------------------------------------------------------------------------
!1. Dependiente del ejemplo (denotado por etiqueta !ul): Dar valores de velocidad y longitud caracter�sticas en la f�rmula del n�mero de Re para calcular 
! Re (formulas escritas en c�digo para cada ejemplo que fue utilizado).
!2. Dependiente del ejemplo (denotado por etiqueta !ac): Dar buena condici�n inicial o aproximaci�n inicial (CI-aprox) para flujo superficial y 
! evitar definir vibv(u)=z(u) para flujo subterr�neo. Dar s�lo valor trivial si navier='si' o navier='no',ap='si' (altura=calado=0). De todos modos
! la selecci�n previa evita que existan calados negativos o nulos los cuales no permitir�an soluci�n o generar�an n�meros de valor infinito (Nan). 
!3. Dependiente del ejemplo: Selecci�n de dominio superficial a trav�s de propiedad, modificaci�n de hmin (abrir o cerrar las etiquetas !es). 
!4. Mala edici�n de la malla: mal dadas las CC en los ficheros mallainicial.txt o mallasubinicial.txt (por ejemplo hay CC de vel=0 entre nodos con CC de 
! H=altura de la l�mina en mallainicial.txt, falta de CC en alg�n nodo del contorno,...). Parece que en el modelo conjunto no es posible prescindir de CC 
!de nivel para el flujo subterr�neo en estacionario (al igual que no lo es el sup aislado estacionario). 
!*** Cambiada la �ltima frase que hab�a por la �ltima frase que hay en el p�rrafo de arriba. Ahora no hay niveles en el contorno m�vil.
!5. Dependiente del ejemplo: faltan ficheros. El fichero mallainicial es necesario para los modelos superficial y conjunto. El fichero mallasubinicial es 
! necesario para los modelos subterr�neo y conjunto. El fichero propiedad.txt es por defecto necesario (a veces no se dan todas las propiedades con �l). 
! En caso de no utilizarlo comentar l�neas (con la etiqueta !fi). Siempre considerar dar propiedades constantes (en las l�neas con la etiqueta !fa),
! propiedades no nulas con sentido f�sico son necesarias para el modelo subterr�neo.
!6. No convergencia: La variable tol1=tol dada en la subrutina gradientesbiconjugados puede representar problemas. En la resoluci�n del sistema lineal no 
! se debe alcanzar un error superior al que se pide en los modelos tol2=tol (del sistema no lineal). Si se alcanza se impidir�a la convergencia y este error 
! suele ser mayor a tol1 (error estimado<tol1 pero error>tol1) en alguna componente. Por ello tol1 debe ser m�s peque�o.     
! Se ha visto que en los peores casos llega con tol1=1e-10 para tol2=1e-6 (caso programado) o tol1=1e-9 para tol2=1e-5. Adem�s en la representaci�n se 
! diferenciar�n mejor las variaciones aunque no est� resuelto el problema no lineal con gran precisi�n.
!7. Evitable en ejecuci�n (usuario elige por pantalla): Se est�n usando las integrales de contorno de t�rminos viscosos (tambi�n al usar el m�todo 
! de Newton) y esto da problemas.
!8. Evitable en ejecuci�n (usuario elige por pantalla): Si no se logra soluci�n para el sistema lineal, el precondicionador para resolverlo puede no 
! ser adecuado y se aconseja elegir el otro (a trav�s de PCGB - subrutina gradientesbiconjugados � PCGBLU - subrutina dslubc, ser� m�s r�pido el LU). 
! En el modelo subterr�neo cualquier precondicionador funcionar� bien. En el superficial: el buen funcionamiento de cada precondicionador depender� del 
! ejemplo. 
! No se da la opci�n del uso del precondicionador LU en caso de usar el m�todo de Newton o la estabilizaci�n (modificar el c�digo para ello) porque se 
! han visto peores resultados.
!9. Improbable. El precondicionador diagonal podr�a no funcionar bien en el modelo superficial porque podr�an aparecer n�meros infinitos (Nan en la 
! variable di). De todos modos, est� programado para que no suceda (ver sub nuevamalla en el caso particular de considerar ciertos elementos).  
!10. No convergencia: En caso de que el modelo superficial no converja (para Re altos pasa) considerar usar estabilizaci�n (sistema mejor condicionado). 
! �sta es adecuada si los elementos de la malla tienen forma regular (misma forma, pudiendo tener diferente tama�o) y funcionar� mejor tanto m�s similar 
! sea la forma.
!11. No convergencia del m�todo de Newton: Si se est� usando el m�todo de newton y no converge tal vez llegue con subir el n�mero de iteraciones previas 
! de Picard en el c�digo. 
!12. Imposible. Al usar el esquema estabilizado en transitorio, los t�rminos de masa con pesos SUPG-PSPG generan problemas, sobre todo si se dan incrementos 
! de tiempo peque�os (los t�rminos van dividos por estos incrementos). Cambiando el c�digo (partes comentadas) se pueden utilizar.
!13. No convergencia inducida: hmin=1e-3 podr�a no ser suficientemente grande como para evitar problemas de localizaci�n y deslocalizaci�n de los mismos 
! elementos indefinidamente. 
!14. Soluciones con error o no convergencia: Si se aplica el modelo superficial con un esquema semi-impl�cito utilizar un incremento de tiempo peque�o o 
! habr� errores en la soluci�n (para valores suficientemente altos no habr� convergencia). Tampoco se recomienda usar un valor muy peque�o porque tardar�a 
! demasiado (a�n no se ha limitado el incremento de tiempo con el par�metro CFL).
!15. Si va a haber (al menos) una zona de flujo superficial aislado es necesario utilizar simulaci�n transitoria. Tampoco se pueden utilizar iteraciones 
! previas de NS 2D en este caso.
!16. Si se utiliza la posibilidad comentada en !c� en mallaaguassomeras se produce no convergencia (comentado all� lo que hacer para usarla).
!17. Si se utilizan iteraciones previas de NS 2D habr� que seleccionar una nu con la que se tenga un gradiente parecido al de la soluci�n.

!Notas sobre la aprox-CI y la selecci�n de dominios:
!---------------------------------------------------
!1. El fichero mallainicial.txt lleva para el contorno de todo el dominio unas CC donde hay flujo superficial y CC nulas donde no lo hay (CC superficiales). 
! El fichero mallasubinicial.txt lleva para el contorno de todo el dominio unas CC donde hay flujo subterr�neo y CC nulas donde no lo hay (CC subterr�neas). 
!2. NS2D 
! La ecuaci�n no tiene en cuenta cota de fondo variable y el proceso iterativo permite calados negativos.
! Est� dada una aproximaci�n inicial con soluci�n trivial para estacionario, le vale cualquiera.
! El transitorio necesita de una CI=vib (uso de vib) y tambi�n vale la soluci�n trivial. 
!3. Aguas someras 
! Poner aproximaci�n para arrancar el estacionario (uso de vib). 
! Con ella, adem�s se efect�a una selecci�n de dominio pues el proceso iterativo no permite calados negativos (luego se hace cargo la condici�n seco-mojado).
! Otra opci�n es usar una soluci�n tras una (o m�s) iteraci�n de la ec NS2D (ap='si'). Aqu� no es necesaria una aprox espec�fica (se usa vib) y se efectuar�
! la selecci�n del dominio del mismo modo.
! En caso donde la l�mina de agua pueda ser pr�cticamente horizontal (CC(H)>zmax de zona seleccionada) vale cualquier opci�n. 
! Es razonable una aprox con altura contante y velocidades nulas que provoque calados. 
! En otro caso, para la primera opci�n se hace necesario tener una soluci�n real (sol con calados y velocidades, lo que es imposible) por lo que se usar� 
! la segunda opci�n con un Re que permita calados en un dominio seleccionado. Se recomienda seleccionar el dominio a trav�s de una propiedad tipo 
! coef de manning (es casi imposible que una aprox-CI inventada permita una selecci�n realista).
! Truco: subir hmin en nuevamalla (c. seco-mojado), hacer m�s it previas de NS (por ejemplo hasta que no haya variaci�n de elementos de dominio superficial 
! entre 2 it), buscar con N-S 2D un calado constante (con NS2D las velocidades son conservativas en 2D indep del calado y para un calado cte lo ser�n en 3D),
! cuidar que no se seque la zona principal con flujo por c. seco-mojado, usar un coef. Manning peque�o (o nulo).  
! Poner una CI para arrancar el transitorio (uso de vib). S�lo se podr� utilizar N-S 2D (ap='si') en el primer incremento (se usa vib).
!4. Aguas someras con uso de aprox inicial con altura constante (casos donde no sea necesaria una selecci�n a trav�s de una propiedad): si se resuelve 
! el caso estacionario con el modelo superficial o el conjunto (con aguas someras), y no se quiere resolver inicialmente en todo el dominio para las 
! CC superficiales se necesita una aprox inicial con menor cota para la altura de agua. 

!Modelo conjunto (para cada incremento de tiempo o para el transitorio se empieza por el modelo superficial):
!------------------------------------------------------------------------------------------------------------
!1. Condiciones de contorno: Se deben diferenciar las CC de flujo superficial de las CC de vel subterr�neas aplicadas como CC superficial.
! Para simular flujo superficial se ve la necesidad de imponer velocidades superficiales y calados. 
! Las citadas velocidades subterr�neas apenas generar�n flujo superficial (no se obtiene el flujo propio de un r�o con CC de calado y 
! condiciones de este tipo, de infiltraci�n). Con ellas es posible simular zonas de agua embalsada (flujo superficial aislado), habiendo 
! CC de vel subterr�nea �nicamente en el contorno m�vil. Adem�s la soluci�n tendr� sentido f�sico (si s�lo hay infiltraci�n).
!2. Flujo aislado: El modelo no genera zonas aisladas superficiales nuevas, s�lo tiene en cuenta las zonas aisladas superficiales localizadas 
! previamente o las localizadas dentro del subdominio superficial (por ejemplo si en una simulaci�n transitoria la altura de agua decrece). 
! El flujo superficial puede generar subdominios subterr�neos (flujo subterr�neo aislado) durante el proceso iterativo.
! El flujo subterr�neo no puede generar subdominios superficiales (flujo superficial aislado) durante el proceso iterativo.
! Las zonas de flujo superficial aislado se localizar�n con la CI-aprox o una propiedad (p.e. n�mero de manning) y si se deslocalizan no se vuelven 
! a tener en cuenta (la consideraci�n de todo el dominio por la CI-aprox obligar� a que las alturas se ajusten a la altura definida como CC superficial).
! Si existe una zona as�, ah� el agua debe encontrarse contenida por una depresi�n del terreno para asegurar que no haya vertido 
! (seleccionar una buena CI-aprox). Si esto no se cumple (puede ocurrir para un tiempo de simulaci�n determinado), la zona aislada se acabar� uniendo 
! a otra zona de forma err�nea ya que la condici�n seco-mojado sumar� elementos sin considerar el incremento de tiempo. Adem�s no se tendr�n condiciones 
! de flujo superficial en el contorno de la zona, el agua no tendr� la suficiente velocidad como para que i sea m�s o menos I existiendo 
! calado normal, por lo que no se obtendr� una soluci�n del tipo de un afluente bajando a un r�o.
!3. Se ha visto que el modelo conjunto no permite el c�lculo de una soluci�n estacionaria si se est� simulando flujo superficial aislado (sin CC de flujo 
! superficial en el contorno de todo el dominio). Tal vez porque la soluci�n estacionaria conlleve al vertido (y sean necesarias esas CC).
!4. La exactitud en la velocidad subterr�nea impuesta en el contorno m�vil depender� del tama�o de los elementos interfaz pegados al contorno ya que
! es calculada a trav�s de un nivel fre�tico medio (un plano). El caudal que se intercambia entre los modelos es funci�n de ella.
!*** Alfa es un coeficiente de calibraci�n en la condici�n de goteo.
 
!Posibilidades (etiquetas !c�):
!------------------------------
!1. Existen dos posibilidades para dar CC (intersecci�n contorno m�vil-contorno) en el modelo conjunto desde la subrutina mallaaguassomeras.
!2. Postproceso existen dos formas de calcular el calado en nodos cuadr�ticos pertenecientes al dominio subterr�neo cercanos al contorno m�vil.
!3. Posibilidad de calcular la soluci�n subterr�nea con una matriz de masa concentrada (te�ricamente mejora el sistema y no tiene error en el 
!caso estacionario).
!4. Con la variable yn se podr�a calcular una soluci�n para aguas someras con ix (pend en x)=Ix e iy (pend en y)=Iy independiente de la velocidad. 
! Siempre generar� soluciones con el mismo calado en todos los nodos independientemente de la cota del terreno con cierto sentido en el caso de 
! canales con igual cota de fondo en cada secci�n transversal (canal rectangular). En este caso se generar�a un calado normal para cualquier velocidad. 
!5. C�lculo preproceso de pendientes medias (a trav�s de los diferentes valores nodales que se tienen) en direcciones x e y multiplicadas por una 
! constante. No utilizada.

!-------------------------------------------------------------------------------------------------------------------------------------------------------
!INICIO DE PROGRAMA
!-------------------------------------------------------------------------------------------------------------------------------------------------------

!M�dulo para los vectores de mayor dimensi�n
!---------------------------------------------
!Los m�dulos permitir�n generar variables globales que pueden ser utilizadas y modificadas desde cualquier subrutina.
module allocatacion
!Al igual que en ciertas subrutinas, se utiliza dimensionamiento din�mico (se dar� la dimensi�n en tiempo de ejecuci�n con allocate).
!Las siguientes ser�n las variables que usan m�s memoria. Se han usado las m�nimas posibles.
integer*4, dimension(:),allocatable::ita,cia,cja 			 
real*8, dimension(:),allocatable::sa,ca
endmodule

!M�dulo para las subrutinas donde se produce intercambio de condiciones de contorno
!------------------------------------------------------------------------------------
module interaccion
integer*4 b,c,e,f,sb(4),eaf,no(6)								  
character a*80
endmodule

!M�dulo para las subrutinas donde se calculan las matrices elementales
!-----------------------------------------------------------------------
module elemental
integer*4, dimension(:),allocatable::pos,posi,posdosi,postresi
integer*4 u,uu,uuu,ui,uj,n(6)								  
real*8 xa,ya,xb,yb,luno(7),ldos(7),ja(7),lunot(7),ldost(7),jat(7),jac
character ac*80
endmodule

!------------------------------------------------------------------------------------------
!Programa principal
!'i' es el n� de nodos, 'j' el n� de elementos, 'modelo' gestiona que modelo se resuelve	 
!------------------------------------------------------------------------------------------
program flujo2D
!Se usa un m�dulo como se har� en otras subrutinas.
use elemental
!Se definen variables globales (programa principal).
integer*4, dimension(:),allocatable::vn 
integer*4 i,j 
character modelo*12	

12  format(/,I5)
23  format(4/,A80)
26  format(6X,6(X,I5))

!Fichero donde se escribir� la soluci�n	y fichero donde se escribir�n las propiedades.
open(unit=8,file='C:\solucionaguassomeras-subterraneo.dat',status='unknown')
open(unit=10,file='C:\solucionpropiedades.dat',status='unknown')	 

!Selecci�n del modelo a utilizar.
!En caso de ser necesaria la condici�n seco-mojado para el modelo superficial o la condici�n similar para el modelo conjunto
!se trabajar� con mallas menores a la malla de todo el dominio (en los ficheros malla.txt y mallasub.txt) para el c�lculo de los coeficientes
!generando s�lo las matrices elementales necesarias. En este caso el programa trabajar� siempre con la misma numeraci�n de nodos y se dar�n CC nulas 
!fuera de la malla para construir un sistema menor.
write(6,*)'Escribe el modelo a usar: superficial/subterraneo/conjunto'
read(5,*) modelo
if (modelo.eq.'superficial') then
 write(8,'(A)')'TITLE = "Aguassomeras"'
 write(8,'(A)')'VARIABLES =  X, Y, Htzt, H, VELX, VELY, VEL, TENX, TENY,'
 write(8,'(A)')'TENXY, VORT' 
 write(6,*)'Debe estar en C:\ el fichero mallainicial.txt'	 
 open(unit=1,file='C:\mallainicial.txt',status='old')
elseif (modelo.eq.'subterraneo') then
 write(8,'(A)')'TITLE = "Subterraneo"'
 write(8,'(A)')'VARIABLES =  X, Y, Hd, E, VELX, VELY, VEL'	  
 write(6,*)'Debe estar en C:\ los ficheros mallasubinicial.txt'	 
 open(unit=1,file='C:\mallasubinicial.txt',status='old')
else
 write(8,'(A)')'TITLE = "Aguassomeras y subterr�neo"'
 write(8,'(A)')'VARIABLES =  X, Y, HtHd, Htzt, H, E, VELX, VELY, VEL' 
 write(6,*)'Deben estar en C:\ los ficheros mallainicial.txt, mallasubinicial.txt'	
 open(unit=1,file='C:\mallainicial.txt',status='old')
endif

!Lectura de n�mero de nodos y elementos de la malla desde el fichero existente.
read (1,12) i
read (1,12) j
rewind(1)

!En los vectores pos se guardar� la posici�n donde empieza la siguiente caja para cada fila.
allocate (vn(i),pos(i),posi(3*i),posdosi(3*i),postresi(3*i))	 

do u=1,i
vn(u)=0
enddo
read(1,23)ac 
 do u=1,j
  !Si vn es 1 el nodo es cuadr�tico. As� vn llevan los nodos "no esquina" de toda la malla.
  read(1,26)n(1),n(4),n(2),n(5),n(3),n(6)	  
  do uu=1,3
  vn(n(3+uu))=1
  enddo
 enddo 
close(1)
	  
call aguassomerassubt (i,j,vn,modelo)
deallocate(vn,pos,posi,posdosi,postresi)
close(8)
close(10) 
end

!-----------------------------------------------------------------------------------------------------------------------------------------------------------
!Subtrutina AGUASSOMERASSUBT donde se calcula la soluci�n para flujo superficial, subterr�neo o conjunto.
!'i' es el n� de nodos, 'j' el n� de elementos, l el n�mero de nodos esquina (de la toda la malla), 
!'ma' el n�mero de Manning, 'nu' la viscosidad (le�das de fichero en 'zzn'), 're' el n�mero de Reynolds, 'x,y' las coordenadas 2D de cada nodo, 
!'z' la cota del terreno en cada nodo, 'zp' la cota del sustrato impermeable en cada nodo, 'zzp' esta cota pero modificada en el contorno m�vil, 
!'kx,ky,ag,nd' las permeabilidades en las direcciones x e y, el �ngulo de anisotrop�a y la porosidad del terreno (le�das de fichero en 'zz,zzz,zza y zzzz'),
!'qx qy' son los caudales para el flujo subterr�neo, 'velx,vely' las velocidades del flujo subterr�neo (obtenidas postproceso),    
!'vn,ve,ven,vo,vnv,vov' son contadores de nodos cuadr�ticos, lineales o con condici�n de contorno ('vnv,vov' para flujo subterr�neo), 
!'vt,vb,vtv,vbv' son los vectores donde se guarda la soluci�n del sistema ('vtv,vbv' para flujo subterr�neo),
!'vit,vib,vitv,vibv' son los vectores soluci�n obtenidos en el instante de tiempo anterior ('vitv,vibv' para flujo subterr�neo),
!'sino' permite la elecci�n del caso estacionario o no estacionario,
!'At,Ata,tac' llevan el valor del incremento de tiempo, el valor de tiempo simulado antes de imponer m�s tiempo y el valor total de tiempo que se simular�,
!'tiempo,nt' el tiempo de simulaci�n y el n�mero de incrementos de tiempo que supone,
!'newton,bcg' permite seleccionar los m�todos de Netwon/Picard, el solver BCG con precondicionador diagonal o LU,
!'navier,est,imp' permite seleccionar las ecuaciones N-S 2D o Aguas someras (modelo superficial), formulaci�n BG o SUPG-PSPG, impl�cito o semi-impl�cito,
!'ap' permite usar la ecuaci�n de N-S 2D en las primeras iteraciones cuando se usan las ecuaciones de aguas someras,
!'ten' permite el c�lculo de las integrales de contorno de los t�rminos viscosos, 'modelo' gestiona que modelo se resuelve, 
!'vib,vit,vibv, vitv' guardan los valores de la iteraci�n anterior en el tiempo (vibv,vitv para flujo subterr�neo),
!'qb' lleva las condiciones de contorno de bombeo para flujo subterr�neo, 
!'ql' lleva las condiciones de contorno de lluvia para flujo subterr�neo y superficial.
!'te,vor' la tensi�n y vorticidad calculada postproceso (s�lo con los resultados del modelo superficial).
!'it,itt' llevan el n� de iteraci�n del m�t. para resolver la no linealidad de las ec. flujo superficial (Newton o Picard) y el n� de la iteraci�n temporal,  
!'u,ui,ut' contador para bucles generales y contador para bucle temporal. 
!'frec,ni' permite la impresi�n de resultados cada cierto n�mero de incrementos de tiempo. 
!-----------------------------------------------------------------------------------------------------------------------------------------------------------
subroutine aguassomerassubt (i,j,vn,modelo)
use interaccion
!Se definen nuevas variables, que ser�n locales (no tra�das desde otra subrutina o el programa principal). Las variables locales deber�an
!desaparecer (de la memoria) al terminar una subrutina (aunque aqu� ya se hace a trav�s del dimensionamiento din�mico para los vectores). 
integer*4 i,j,u,ui,ut,vn(i),it,itt,ss,l,nt,ni,frec,app				 		  	 
real*8, dimension(:),allocatable::x,y,z,zp,zzp,vt,vit,vib,vb,vo,eval,vov,vtv,vitv,vibv
real*8, dimension(:),allocatable::kix,kiy,ag,nd,ma,qx,qy,qb,ql,velx,vely,vx,vy,ht,evol,te,vor  
real*8 zu,tol,re,nu,mu,At,Ata,tiempo,tac,zz,zzz,zza,zzzz,zzn,qbb,qll  
character modelo*12,sino*2,newton*2,bcg*2,ten*2,Atee*10,Atai*29,ww*1,navier*2,est*2,imp*2,ap*2
!Se utiliza dimensionamiento din�mico aunque no sea necesario (se podr�a definir como la variable tra�da vn, que adem�s es global). 
!Sin embargo es el modo correcto de proceder para para que la subrutina pre-asigne la memoria necesaria para las nuevas variables locales.
allocate(x(i),y(i),z(i),zp(i),zzp(i),vt(3*i),vit(3*i),vib(3*i),vb(3*i),vo(3*i),eval(i),vov(i),vtv(i),vitv(i),vibv(i))
allocate(kix(i),kiy(i),ag(i),nd(i),ma(i),qx(i),qy(i),qb(i),ql(i),velx(i),vely(i),vx(i),vy(i),ht(i),evol(i),te(3*i),vor(i))
	  			 											   
!Formatos de lectura o escritura:
!--------------------------------
15  format(6X,3(2X,F12.4))			 
17  format(X,I5,5(4X,F14.12))
18  format(X,I5,X,A1,X,F11.7)		  
21  format(3/,I5)
!Se escribir�n el n� de elementos equivalentes (habr� un m�ximo de 399996 elementos P2-P1, por lo que llega con 6 cifras), de tres nodos.
22  format('ZONE T=',A29) 
23  format('I=',4x,I5,',',2x,'J=',4x,I6,',',x,'    DATAPACKING = POINT, ZONETYPE = FETRIANGLE')	  
24  format(2x,'ZONE  I=',4x,I5,',',2x,'J=',4x,I6,',',x,'    DATAPACKING = POINT, ZONETYPE = FETRIANGLE')
!As�, se escribir�n tres columnas con a b c (siendo a,b,c los tres nodos de cada elemento).
!Tambi�n vale: ...,x,'    F=FEPOINT') pero luego hay que escribir cuatro columnas con a b c c cambiando el formato 41 a format(4(3x,I5)).
!Necesario representarla con FEPOINT si se representa en Tecplot 9.0. En versiones siguientes es mejor hacerlo con FETRIANGLE (ahorra espacio de fichero). 																											  
26  format(6X,6(X,I5))
27  format(4/,A80)
!Para el fichero de resultados:
!Se escribir�n las coordenadas con el mismo formato con el que se leen, con 7 enteros (o 6 y signo negativo) y 4 decimales.
!Se escribir�n las propiedades con 4 enteros(o 3 y signo negativo) y 10 decimales.
38  format(2(x,F12.4),5(2x,F15.10))  
39  format(2(x,F12.4),5(2x,F15.10))   
40  format(2(x,F12.4),7(2x,F15.10))  
41  format(3(3x,I5))
48  format(3(x,F12.4),(2x,F15.10))	  
49  format(3(x,F12.4),3(2x,F15.10))  
50  format(4(x,F12.4),4(2x,F15.10))  
	  
!Inicializaci�n de variables:
!----------------------------
!Se utiliza doble precisi�n para los c�lculos y es necesario escribir los n�meros reales de forma que se considere
!doble precisi�n. As�, se utiliza la notaci�n 0.1125_8, que es equivalente a 1.125*1d-1, para escribir el n�mero 0.11250000000...		
itt=0
Ata=0.0_8
tac=0.0_8
do u=1,3*i
vt(u)=0.0_8
vb(u)=0.0_8
enddo
do u=1,i
z(u)=0.0_8
zp(u)=0.0_8
zzp(u)=0.0_8 
vtv(u)=0.1_8
velx(u)=0.0_8
vely(u)=0.0_8
eval(u)=0.0_8
ql(u)=0.0_8
qb(u)=0.0_8
!Valores para la primera iteraci�n del m�todo para resolver la no linealidad.
!Los siguientes valores s�lo se usan en caso de usar el modelo superficial con N-S 2D (valen valores nulos) o el modelo subterr�neo (no valen).
!En este caso, al buscar sol. estacionaria o transitoria se usa la solucion trivial (vt=0) para N-S 2D. 
!En otro caso, al buscar sol. estacionaria o transitoria, se escribe despu�s el valor de la aproximaci�n inicial 
!o la CI (estar� en vib), tal que vit=vib-z/vit=vib,vb=vib,vitv=vibv.
!Valores de velocidad para flujo superficial.
vit(u)=0.0_8
vit(i+u)=0.0_8
!Valores de espesor (flujo superficial) y altura (flujo subterr�neo).
vit(2*i+u)=0.0_8			   
vitv(u)=0.1_8
enddo
!Se podr�a usar tambi�n el comando data para 'sb' con el que se definen datos que no cambian. No es el caso de b,c,e y por tanto 
!no se podr�a definir data b/0/,c/0/,e/0/ (equivalente a: data b,c,e/0,0,0/). 
sb=(/0,0,3,3/)
b=0
c=0
e=0
f=0
nt=0
ni=1

!Selecci�n de estabilizaci�n, tipo de esquema y ecuaciones superficiales a utilizar.
!-----------------------------------------------------------------------------------
!Si se resuelve con el modelo subterr�neo siempre se utiliza un esquema impl�cito.
if (modelo.eq.'superficial')then
 write(6,*)'Modelo Navier Stokes en vez de aguas someras? (si/no)'
 read(5,*)navier
 !Resoluci�n con/sin estabilizaci�n para aguas someras y para Navier-Stokes 2D.
 write(6,*)'Resolucion con esquema estabilizado en vez de Bubnov-Galerkin? (si/no)'
 read(5,*) est
 if ((est.eq.'si').and.(navier.eq.'no')) then
 write(6,*)'Ojo, no se aplica lluvia (no programado)'
 endif
 !Resoluci�n con esquema impl�cito/semi-implicito (incondicionalmente estables) para aguas someras y para Navier-Stokes 2D.
 !C�lculo de soluci�n estacionaria: se llegar� a las mismas soluciones si CC constantes	y es posible una soluci�n estacionaria.
 !Con impl�cito (resoluci�n transitoria) o semi-impl�cito se puede calcular con incrementos de tiempo (tambi�n si fuese expl�cito).
 !En impl�cito tambi�n se puede calcular estacionario ignorando t�rminos temporales (resoluci�n estacionaria).   
 !C�lculo de soluci�n transitoria: En ambos casos se requiere soluci�n con incrementos de tiempo.
 write(6,*)'Resolucion con esquema implicito en vez de semi-implicito? (si/no)'
 read(5,*) imp
elseif (modelo.eq.'conjunto')then
 !Siempre se usar� aguas someras.
 navier='no'
 !Resoluci�n con/sin estabilizaci�n para aguas someras.
 write(6,*)'Resolucion con esquema estabilizado en vez de Bubnov-Galerkin? (si/no)'
 read(5,*) est
 if (est.eq.'si') then
 write(6,*)'Ojo, no se aplica lluvia en modelo aguas someras (no programado)'
 endif
 !Siempre se usar� un esquema impl�cito.
 imp='si'
endif 
ap='no'
app=0

!Lectura de propiedades:
!------------------------
!C�lculo de la variable l
l=i-count(vn.eq.1)
!fi (hasta close(12) inclusive)
eaf=0
open(unit=12,file='C:\propiedad.txt',status='old',iostat=eaf)
if (eaf.eq.0)then
read(12,'(A)')a    
 !Valores de conductividad en m/s, de �ngulo de anisotrop�a en radianes, de porosidad (adimensional) y de manning.
 !Lectura de propiedades desde fichero (conductividad y porosidad necesarias para flujo subterr�neo) 
 do u=1,l	   
  read(12,17)ss,zz,zzz,zza,zzzz,zzn
  kix(ss)=zz
  kiy(ss)=zzz
  ag(ss)=zza
  nd(ss)=zzzz
  ma(ss)=zzn
 enddo
else				 
write(6,*)'Ojo, falta propiedad.txt modificar codigo segun sus especificaciones de error'
endif
close(12)

!Lectura de coordenadas y condiciones de contorno para el modelo superficial:
!----------------------------------------------------------------------------
!L�nea a l�nea del fichero, empieza a leer las coordenadas x,y,z de cada nodo 
!en la l�nea n�elementos+6 del fichero para flujo superficial.
if ((modelo.eq.'superficial').or.(modelo.eq.'conjunto')) then 
open(unit=1,file='C:\mallainicial.txt',status='old')
do u=1,6+j
read(1,'(A)')a
enddo
 !Valores de coordenadas en metros (incluyendo la cota del terreno).
 do u=1,i
 read(1,15) x(u),y(u),z(u) 
 !fa (hasta ma(u) inclusive)	
 !Si no existe fichero de propiedades no se da propiedad alguna, y se debe dar un valor para Manning (vale dar valor nulo) para inicializar la variable.
 !Si existe puede no haberse dado esta propiedad (si todos los valores le�dos son nulos).
 !Si existe y se tienen valores para la propiedad aqu� se pueden sobreescribir. 
 !Un valor constante con sentido f�sico (para hormig�n):
 ma(u)=0.014_8 
 enddo
read(1,'(A)')a
do u=1,3*i
vo(u)=sqrt(2.0_8)
enddo
!Lectura de las condiciones de contorno para flujo superficial desde fichero.
!El programa lee ficheros donde hay "i" l�neas con condici�n para la velocidad en x, e "i" l�neas en y. Se permite cualquier n� l�neas con condici�n 
!para el calado (por ejemplo tantas como nodos esquina como se tiene por defecto). 
!Se obtiene el vector de valores conocidos para cada iteraci�n (CC Dirichlet constantes).
u=0
eaf=0
do while (eaf.eq.0)
 u=u+1	   
 read(1,18,iostat=eaf)ss,ww,zu
 if (eaf.eq.0) then			
  if (ww.eq.'o') then
    cycle
  elseif (u.le.i) then
    vo(ss)=zu
  elseif ((u.gt.i).and.(u.le.(2*i))) then
    vo(i+ss)=zu
  else
    vo(2*i+ss)=zu+z(ss)
  endif
 endif 
enddo
!Se guarda el fichero con toda la malla 'mallainicial.txt' en otro archivo (para no modificar los ficheros de datos originales si se generan
!mallas de menos elementos). El nuevo archivo ser� el que se utilize durante la simulaci�n.
rewind(1)
open(unit=4,file='C:\malla.txt',status='unknown')
 do u=1,5+j
 read(1,'(A)')a 
 write(4,'(A)')a  
 enddo
close(4)
close(1)
!Valor del n�mero de Reynolds dado por pantalla. 
write(6,*)'Introduce un valor para el numero de Reynolds: '
read(5,*) re
!Calculo de la viscosidad para ese n�mero de Reynolds (nu=vel*L/re). Depende del ejemplo.
!'vel' es la velocidad caracter�stica (para calcular la viscosidad). Ser� el valor de la velocidad promedio condicion de contorno. 
!'L' es la longitud caracter�stica. Ser� el ancho en la direcci�n perpendicular a la direcci�n de avance del fluido. 
!En ejemplos utilizados: nu=(1.28_8*30.0_8)/re  nu=(0.256_8*0.43_8)/re  nu=(1.0_8*400.0_8)/re  nu=(0.666666*2.0_8)/re  nu=(0.33_8*6.0_8)/re   !ul
nu=(1.28_8*30.0_8)/re	
mu=nu
 !Resoluci�n con Newton/Picard para las ecuaciones de aguas someras con esquema impl�cito y formulaci�n BG. 
 if ((est.eq.'no').and.(navier.eq.'no').and.(imp.eq.'si').and.(modelo.ne.'subterraneo')) then
 write(6,*)'Quieres aplicar el metodo de newton para el modelo superficial?(si/no)'
 read(5,*)newton
 else
 newton='no'
 endif
endif

!Lectura de coordenadas y condiciones de contorno para el modelo subterr�neo:
!----------------------------------------------------------------------------
!L�nea a l�nea del fichero, empieza a leer las coordenadas x,y,zp de cada nodo 
!en la l�nea n�elementos+6 del fichero para flujo subuterr�neo.
if ((modelo.eq.'subterraneo').or.(modelo.eq.'conjunto')) then         
open(unit=3,file='C:\mallasubinicial.txt',status='old')
do u=1,6+j
read(3,'(A)')a
enddo
 !Valores de coordenadas en metros (incluyendo cota del sustrato)
 do u=1,i
  read(3,15) x(u),y(u),zp(u) 
  !fa  (hasta nd(ss) inclusive)
  !Si no existe fichero de propiedades no se da propiedad alguna, y se deben dar valores de kix,kiy,ag y nd para inicializar las variables.
  !Si existe puede no haberse dado alguna propiedad (todos los valores le�dos son nulos).
  !Si existe y se tiene una propiedad aqu� se pueden sobreescribir. 
  !Propiedades de valor constante con sentido f�sico (p.e.):
  !kix(ss)=0.00001_8
  !kiy(ss)=0.00001_8
  !ag(ss)=0.0_8
  !nd(ss)=0.01_8
 enddo 
read(3,'(A)')a
do u=1,i
vov(u)=sqrt(2.0_8)
qx(u)=sqrt(3.0_8)					 
qy(u)=sqrt(3.0_8)
enddo
u=0
eaf=0
!Lectura de las condiciones de contorno desde fichero.
!En el fichero han de aparecer "i" l�neas para las CC de caudal por metro lineal aunque s�lo se tendr�n valores para los nodos esquina.
do while (eaf.eq.0)
u=u+1
 read(3,18,iostat=eaf)ss,ww,zz
 if (eaf.eq.0) then			   
  if (ww.eq.'o') then
  cycle 
  elseif (u.le.i) then
    qx(ss)=zz
  elseif ((u.gt.i).and.(u.le.2*i)) then
    qy(ss)=zz
  else
    vov(ss)=zz
  endif
 endif 
enddo 						  
!Igual que antes, se guarda el fichero con toda la malla 'mallasubinicial.txt' en otro archivo.
rewind(3)
open(unit=7,file='C:\mallasub.txt',status='unknown')
 do u=1,5+j
 read(3,'(A)')a
 write(7,'(A)')a
 enddo
close(7)
close(3)
!Se imponen condiciones de bombeo desde pantalla (constantes en el tiempo).
write(6,*)'Existe algun bombeo?(si/no)'
read(5,*)sino
if (sino.eq.'si') then
write(6,*)'Introduce nodo y valor de caudal:'
read(5,*) u,qbb 
qb(u)=qbb
endif
endif

!Selecci�n del precondicionador, y del uso de las tensiones viscosas:
!-------------------------------------------------------------------- 
if ((newton.eq.'no').and.(est.eq.'no')) then
!S�lo en los casos donde funciona bien se permite precondicionador LU. As�, s�lo si se usa formulaci�n BG y Picard se entrar�a, lo que 
!siempre ocurrir� para el modelo subterr�neo (aunque est no tiene valor en este caso).
write(6,*) 'Utilizar PBCGLU en vez de PCGB (si/no)?'
read(5,*)bcg
else
bcg='no'
endif
if ((modelo.eq.'conjunto').or.(modelo.eq.'superficial')) then
write(6,*) 'Utilizar tensiones viscosas (modelo superficial) (si/no)?'
read(5,*)ten
endif 

!Condiciones iniciales o aproximaci�n inicial   !ac
!--------------------------------------------------
!En caso de usar el modelo superficial con N-S 2D o el modelo subterr�neo, s�lo ser� necesario definir vib � vibv como CI si se simula en transitorio 
!(se usar�n valores diferentes de CI y de valor para arrancar el m�todo de resoluci�n del sistema no-lineal, a no ser que se escriban los mismos). 
!En este caso habr� aprox inicial v�lida en vit para el estacionario (para arrancar el m�todo). Es as� porque no se har� vit=vib, vitv=vibv. 
!En otro caso, siempre ser� necesario escribir aqu� el valor de vib � vibv. En transitorio ser� la CI (mismo valor de CI y de valor para arrancar el 
!m�todo de resoluci�n del sistema no-lineal). En estacionario ser� la aproximaci�n inicial (para arrancar el m�todo). Es as� porque se 
!har� vit=vib, vitv=vibv. 
!Flujo superficial: la iteraci�n no lineal requerir� soluci�n con calados (vit) y la iteraci�n temporal requerir� soluci�n con alturas (vib).
do u=1,i		        
!Valores de velocidad para flujo superficial.
vib(u)=0.0_8           
vib(i+u)=0.0_8         
!Valores de altura de la l�mina de agua para flujo superficial (en vib) y de altura de agua o nivel fre�tico para flujo subterr�neo (vibv).
!Se deber�a dar el mismo valor inicial para ambos modelos de forma que vib(2*i+u)=vibv(u). 
!Ojo con definir tambi�n vibv si se usa el modelo superficial de aguas someras por la posterior selecci�n de dominios (es recomendable 
!definir vibv=vib si se desea utilizar todo el dominio).
!Ojo si vibv(u)=z(u) porque con la posterior selecci�n de dominios aquellas zonas que no pertenezcan cumplir�n vibv(u)=z(u)+vt(2*i+u) en 
!la subrutina mallaaguassomeras.
!Si no tiene sentido una CI tipo plano con altura constante (lo tiene si z(CCv)>CCH) ser� dif�cil darla. Para la posterior selecci�n del dominio 
!se puede usar Manning y para el c�lculo estacionario no har� falta definir nada si ap='si'. 
!vibv(u)=z(u)-0.1_8 
!if (z(u).le.33.48_8)then
!vib(2*i+u)=34.48_8  	
!else
!vib(2*i+u)=z(u)+1.48_8
!endif
!En ejemplos utilizados (planos con altura constante):
!Un plano de altura constante (tanto para altura como para nivel fre�tico). 										 
!vib(2*i+u)=34.48_8
!vibv(u)=34.48_8
!Un plano de altura constante si z<altura y una superficie con altura mayor que ese valor en otro caso (tanto menor que la del terreno
!tanto mayor su cota). Visto que funciona mejor que un plano con altura constante si las CC de caudal subterr�neo son nulas y se aplica lluvia.
!50.0_8 es un par�metro de ajuste de la superficie.
!if (z(u).lt.34.48_8)then	 
vib(2*i+u)=11.50_8
vibv(u)=11.50_8 	
!else
!vib(2*i+u)=z(u)-(z(u)-34.48_8)/50.0_8	 
!vibv(u)=z(u)-(z(u)-34.48_8)/50.0_8
!endif
enddo					

!Selecci�n de subdominios superficiales bien desde la CI o la aproximaci�n inicial (caso estacionario), bien desde el n�mero de manning. 
!---------------------------------------------------------------------------------------------------------------------------------------  
!En caso de usar el modelo superficial con N-S 2D o el modelo subterr�neo, no se realizar� esta selecci�n.
if ((modelo.eq.'conjunto').or.((modelo.eq.'superficial').and.(navier.eq.'no'))) then
do u=1,i
vb(u)=vib(u)
vb(i+u)=vib(i+u)
vit(u)=vib(u) 
vit(i+u)=vib(i+u)
vitv(u)=0.0_8
eval(u)=sqrt(3.0_8) 
!Selecci�n de dominio inicial con la CI o aproximaci�n inicial.
if (vibv(u).lt.z(u)) then
!Selecci�n de dominio inicial con el n�mero de manning (selecc como superf zona con manning 0.05).	 !es
!if (ma(u).ne.0.05_8) then  								                                             !es
evol(u)=0.0_8
else
evol(u)=sqrt(3.0_8)
endif
enddo
!Sobreescritura del coeficiente de Manning tras seleccionar el dominio con �l (hasta enddo). !es
!do u=1,i		 
!ma(u)=0.0_8		
!enddo			
!Obtenci�n de una aproximaci�n a partir de NS2D con selecci�n del dominio posterior a su resoluci�n	(bueno si velocidades nulas en aprox inicial o CI)
write(6,*)'Quieres obtener una primera aproximacion con NS 2D (si/no)'
write(6,*)'Para ello es necesario ajustar la viscosidad'
read(5,*) ap
if (ap.eq.'si') then
write(6,*)'Escribe el numero maximo de iteraciones previas de NS 2D'
read(5,*) app
endif

!Selecci�n del subdominio 																		  
open(unit=4,file='C:\mallainicial.txt',status='old') 
read(4,27)a	 
 do u=1,j
 read(4,26)no(1),no(4),no(2),no(5),no(3),no(6)   
   if ((evol(no(1)).gt.0.0_8).and.(evol(no(2)).gt.0.0_8).and.(evol(no(3)).gt.0.0_8)) then
   do ui=1,3
   vit(2*i+no(ui))=vib(2*i+no(ui))-z(no(ui))
   vb(2*i+no(ui))=vib(2*i+no(ui))
   enddo
   else  
   do ui=1,3
   vitv(no(ui))=vibv(no(ui)) 
   eval(no(ui))=0.0_8
   eval(no(3+ui))=0.0_8
   enddo
   endif 
 enddo
close(4) 
!*** S�lo se dan valores iguales en la primera iteraci�n
do u=1,i						 !*** A�adida esta l�nea
if (vb(2*i+u).eq.vitv(u))then	 !*** A�adida esta l�nea
zzp(u)=z(u)						 !*** A�adida esta l�nea
endif							 !*** A�adida esta l�nea
enddo							 !*** A�adida esta l�nea
!Se observa si existe flujo superficial en todo el dominio o si no existe flujo superficial en todo el dominio.
!Si se resuelve el modelo conjunto se resolver� (inicialmente) s�lo la ecuaci�n de aguas someras o s�lo la ecuaci�n subterr�nea.
if (minval(eval).eq.sqrt(3.0_8)) then 
write(6,*)'...Solo flujo superficial para esa condicion inicial o aproximacion'
elseif (maxval(eval).eq.0.0_8) then
write(6,*)'...No hay flujo superficial para esa condicion inicial o aproximacion'		
endif
write(6,*)' '
endif

!Elecci�n de modelizaci�n estacionaria (esquema impl�cito) o transitoria (esquema impl�cito o semi-impl�cito)
!------------------------------------------------------------------------------------------------------------
!Se construir�n las cajas de masa para las ecuaciones de aguas someras y para la ecuaci�n subterr�nea para el caso transitorio.
!Tambi�n se impone condici�n de lluvia por pantalla, que es considerada para el modelo superficial (si aguas someras y formulaci�n BG) y subterr�neo.
if (imp.eq.'no')then
sino='si'
else
write(6,*)'Quieres obtener la solucion en tiempos en vez de'
write(6,*)'obtener la solucion estacionaria? (si/no)'
read(5,*)sino
endif
do while (nt.eq.0)
 if (sino.eq.'si') then	
 !Simulaci�n transitoria
  write(6,*)' '
  write(6,*)'t=',Ata
  write(6,*)'Escribe el intervalo de tiempo total de simulacion en segundos:'
  read(5,*)tiempo
  write(6,*)'Escribe el incremento de tiempo de resolucion en segundos:'	   
  read(5,*)At
  !Se permite la impresi�n cada cierto n�mero de incrementos de simulaci�n para evitar generar archivos de mucho tama�o (por ejemplo, imprimir cada 
  !2 permitir� generar 5 archivos en vez de 10 al simular tiempo=10 d�as con At=1 d�a). As�, habr� incrementos de tiempo en que no se imprima.
  write(6,*)'Cada cuantos incrementos se imprime?'
  read(5,*)frec
  tac=tiempo+tac
  nt=idint((tac-Ata)/At)
  if (((navier.eq.'no').and.(modelo.ne.'subterraneo')).or.(modelo.eq.'subterraneo')) then
  !El volumen de agua por lluvia ser� multiplicado por At (para cada modelo).
  write(6,*)'Introduce un valor de intensidad media de lluvia'
  write(6,*)'durante el intervalo de simulacion en mm/d (habitual 0-100)'
  endif
 else 
 !Simulaci�n estacionaria
  nt=1
  if (((navier.eq.'no').and.(modelo.ne.'subterraneo')).or.(modelo.eq.'subterraneo')) then
  !El volumen de agua por lluvia ser�a infinito para cualquier intensidad.  
  write(6,*)'Introduce un valor de intensidad de lluvia en mm/d'
  write(6,*)'Debe ser constante en el tiempo (habitual 0-100)'
  endif
 endif
 if	(((navier.eq.'no').and.(modelo.ne.'subterraneo')).or.(modelo.eq.'subterraneo')) then
 read(5,*) qll
 do u=1,i
 !Intensidad de lluvia en m/s (1mm=1L/m2)
 ql(u)=qll/(86400.0_8*1000.0_8)	
 enddo
 endif
do ut=1,nt
!Problema no lineal superficial:
!En aguas someras no es posible tener valores nulos o negativos de calado como valores de la iteraci�n anterior. Si aparecen despu�s la condici�n 
!seco-mojado limita el dominio de resoluci�n.
!Para transitorio se utiliza la soluci�n en instante anterior en cada paso de tiempo para comenzar Picard o Newton.
!Se utiliza la condici�n inicial en el primer paso si aguas someras en modelo superficial o conjunto (sino la trivial). 
!Para estacionario se utiliza la aproximaci�n inicial si aguas someras en modelo superficial o conjunto (sino la trivial).
do u=1,3*i
vt(u)=vit(u)     
enddo
!Problema no lineal subterr�neo:
!En ec. subterr�nea tampoco es posible tener valores nulos o negativos de calado como valores de la iteraci�n anterior. Si aparecen despu�s, la condici�n 
!de espesor m�nimo permite calcular coeficientes que no dan problema.
!Para transitorio se utiliza la soluci�n en instante anterior en cada paso de tiempo para comenzar Picard.
!Se utiliza la condici�n inicial en el primer paso en modelo conjunto (sino la soluci�n 0.1).
!Para estacionario se utiliza la aproximaci�n inicial en modelo conjunto (sino la soluci�n 0.1).
do u=1,i
vtv(u)=vitv(u)    
enddo
													   	  
!Gesti�n de las ecuaciones de N-S 2D, aguas someras y agua subterr�nea.
!----------------------------------------------------------------------
!Resoluci�n de los modelos superficial, subterr�neo y conjunto.  
!Respecto a 'vn' y 'vnv', en caso de tratarse del modelo conjunto (vn,vnv) o del modelo superficial (vn) y existe contorno m�vil (aplicaci�n de condici�n 
!seco-mojado o de la condici�n similar para el modelo conjunto), pueden ser diferentes en cada iteraci�n no-lineal de aguas someras. Cada uno llevar� los 
!nodos cuadr�ticos que tienen la malla superficial y la subterr�nea respectivamente.
!En las ecuaciones de aguas someras con Manning=0 y Re1 se obtendr�n soluciones similares a si se usa Manning=n y Re2 con Re2>Re1.

!Modelo superficial.
if (modelo.eq.'superficial') then
 tol=1e-5
 it=0
 !La condici�n seco mojado se aplica dentro de la subrutina aguassomeras: a partir de la segunda iteraci�n 
 !no-lineal se utiliza el resultado de flujo superficial de la iteraci�n anterior con velocidades nulas en el contorno m�vil.
 !Permite a�adir elementos al dominio superficial o desechar elementos del dominio superficial. 
 if ((navier.eq.'no').and.(itt.eq.0)) then 
 !Necesario usar la subrutina mallaaguassomeras si se selecciona dominio previamente con la CI o la aproximaci�n inicial (aqu� tambi�n se selecciona
 !la malla para flujo superficial, algo que despu�s hace la condici�n seco-mojado) 
 call mallaaguassomeras (i,j,z,zzp,vtv,velx,vely,eval,vo,vt)  	 !*** A�adido zzp
  !Aproximaci�n inicial con las ecuaciones N-S 2D. Se puede usar siempre, cuando se ha seleccionado una resoluci�n impl�cita o semi-impl�cita
  !dici�ndose en el segundo caso el n�mero de iteraciones temporales (ambas ecuaciones se resolver�n con el mismo esquema como debe hacerse). 
  if (ap.eq.'si') then                                                                                              
  !Podr�a usarse un dominio menor al de toda la malla al haberse hecho una selecci�n previa (buscar un Re que permita calados en el dominio 
  !seleccionado.
  nu=0.5_8  
  !Estas iteraciones se hacen con Picard ya que N-S 2D s�lo est� programado para Picard.
  !Se har� una selecci�n tras cada resoluci�n con la condici�n seco-mojado.  
  do u=1,app                                                                                                                       
  call aguassomeras(modelo,'si',ap,tol,i,j,x,y,z,it,itt,nu,ma,sino,bcg,ten,'no',At,vn,vo,eval,ql,vib,vt,vb,est,imp,vit)           
  if (tol.ne.1e-5) then
   !No hecho todo el n�mero previo de iteraciones NS si existe convergencia antes.
   !Reescrita la tolerancia tras iteraciones previas si hay convergencia. As�, se empieza por Sw despu�s.
   tol=1e-5
   exit
   endif
  enddo                                                                                                                          
  nu=mu
  endif
 endif 
 !En aguas someras se aplica el m�todo de Picard y el de Newton y se har� una selecci�n tras cada resoluci�n con la condici�n seco-mojado.
 !Esquema semi-impl�cito - iteraci�n para un incremento de tiempo.
 if (imp.eq.'no')then
 call aguassomeras(modelo,navier,ap,tol,i,j,x,y,z,it,itt,nu,ma,sino,bcg,ten,newton,At,vn,vo,eval,ql,vib,vt,vb,est,imp,vit)
 write(6,*)'t=',Ata+At
 else
 !Esquema impl�cito - inicio de iteraciones para un incremento de tiempo o para el problema estacionario.
 do while (tol.eq.1e-5)
 call aguassomeras(modelo,navier,ap,tol,i,j,x,y,z,it,itt,nu,ma,sino,bcg,ten,newton,At,vn,vo,eval,ql,vib,vt,vb,est,imp,vit)  
 enddo
 write(6,*)'Convergencia modelo aguassomeras/navierstokes en iteracion:',it
  if (sino.eq.'si')then
  write(6,*)'At=',Ata+At
  endif
 endif
 !En ejemplos anteriores:
 !Resultados velocidades en linea central del cuadrado para Cavity flow
 !open(unit=11,file='C:\soluciones linea central.txt',status='unknown')
 !228 format (2x,F10.5,2x,F10.5)
 !do u=1,i
 !if (x(u).eq.200) then
 !write(11,228) vt(u),y(u)/400
 !endif
 !enddo
 !close(11)
 call vorticidad (i,j,x,y,nu,vt,te,vor)

!Modelo subterr�neo.
elseif (modelo.eq.'subterraneo') then 
 !Se usa toda la malla.
 !No se aplica ninguna condici�n que permita seleccionar elementos o desechar elementos del dominio subterr�neo.
 !Se aplica el m�todo de Picard.
 !Esquema impl�cito - Dentro de la subrutina aguassubterranea se hacen iteraciones para un incremento de tiempo o para el problema estacionario.
 call aguassubterranea(i,j,x,y,z,zp,zzp,nd,kix,kiy,ag,sino,bcg,At,vn,vov,eval,vb,ql,qb,qx,qy,vibv,vtv,velx,vely)  !*** A�adidos z,zzp,vb
 if (sino.eq.'si')then
 write(6,*)'At=',Ata+At
 endif

!Modelo conjunto. 
!Se decide si existe flujo superficial y subterr�neo, y si existen ambos se aplican alternativamente las ecuaciones de 
!aguas someras y la ecuaci�n de flujo subterr�neo. En este caso se realiza una iteraci�n de las primeras y se resuelve la segunda (en cada 
!iteraci�n del modelo conjunto). Se realizan estas iteraciones hasta convergencia del modelo superficial, y se hace para 
!cada incremento de tiempo (partiendo de la misma condici�n inicial) si se busca una soluci�n transitoria.
!Se resuelve primero el superficial aunque no haya condiciones de contorno de flujo superficial lo que puede ocurrir si se ha seleccionado previamente
!con la CI-aprox una zona de agua embalsada (en cuyo caso ya no se habr�n definido condiciones de este tipo en el contorno de todo el dominio, y como
!no se est�n usando CC variables en el tiempo siempre se usar�n siempre �stas en caso de simular un transitorio...).
elseif (modelo.eq.'conjunto') then   
 tol=1e-6
 it=0
 !La condici�n similar a la condici�n seco mojado se aplica dentro de la subrutina aguassomeras: a partir de la segunda iteraci�n 
 !no-lineal se utiliza el resultado de flujo superficial de la iteraci�n anterior con velocidades subterr�neas en el contorno m�vil.
 !Permite a�adir elementos al dominio superficial o desechar elementos del dominio superficial (y generar flujo subterr�neo aislado).
 !No se aplica ninguna condici�n en la ecuaci�n subterr�nea que permita seleccionar elementos o desechar elementos del dominio 
 !subterr�neo (ni generar flujo superficial aislado), modific�ndose �ste con la modificaci�n del dominio superficial.
 !Inicio de iteraciones conjuntas: 
  do while (tol.eq.1e-6)
  if (maxval(eval).eq.0.0) then 
  !Toda la soluci�n previa es subterr�nea (de acuerdo a la selecci�n previa � a la �ltima soluci�n), por lo que se resuelve impl�citamente el 
  !subterr�neo hasta convergencia (Picard) y se tiene la soluci�n conjunta para este incremento de tiempo, o la soluci�n conjunta estacionaria.	       
   call aguassubterranea(i,j,x,y,z,zp,zzp,nd,kix,kiy,ag,sino,bcg,At,vn,vov,eval,vb,ql,qb,qx,qy,vibv,vtv,velx,vely)	   !*** A�adidos z,zp,vb
   !Las siguientes dos l�neas son equivalentes a escribir: exit
   tol=1.0
   cycle   
  else
  !La soluci�n previa es toda superficial o superficial y subterr�nea. 	 
  !Siempre necesario usar la subrutina mallaaguassomeras previamente para obtener las velocidades subterr�neas a aplicar en el contorno m�vil (aqu�
  !tambi�n se selecciona la malla para flujo superficial algo que ahora no hace la condici�n similar a la condici�n seco-mojado).
  call mallaaguassomeras (i,j,z,zzp,vtv,velx,vely,eval,vo,vt)		 !*** A�adido zzp
   if (minval(eval).eq.0.0) then
   !Hay soluci�n superficial en al menos un elemento (si sol. superficial s�lo en un nodos aislados no se calcula con la ec. aguas someras).
   !Aproximaci�n inicial con las ecuaciones N-S 2D (ambas ecuaciones superficiales se resolver�n impl�citamente).
   if ((it.eq.0).and.(itt.eq.0).and.(ap.eq.'si'))then
   !Podr�a usarse un dominio menor al de toda la malla al haberse hecho una selecci�n previa (buscar un Re que permita calados en el dominio 
   !seleccionado.                                                                                                  
   nu=0.5_8	    
   !Estas iteraciones se hacen con Picard ya que N-S 2D s�lo est� programado para Picard.
   !Se har� una selecci�n tras cada resoluci�n con la condici�n seco-mojado (no la similar a la condici�n seco-mojado).
   do u=1,app                                                                                                                       
   call aguassomeras('superficial ','si',ap,tol,i,j,x,y,z,it,itt,nu,ma,sino,bcg,ten,'no',At,vn,vo,eval,ql,vib,vt,vb,est,imp,vit)   
   if (tol.ne.1e-6) then
	!En este caso es necesario reescribir la tolerancia para que se eval�e bien la convergencia del sistema no lineal de aguas someras y evitar 
	!que se obligue a hacer una iteraci�n conjunta solo
    tol=1e-6
    exit
    endif
   enddo                                                                                                                          
   nu=mu
   endif
   !Se aplica impl�citamente una iteraci�n de la ecuaci�n de aguas someras (Picard o Newton) para este incremento de tiempo, o para el c�culo de la 
   !soluci�n conjunta estacionaria.	Se har� una selecci�n tras la resoluci�n con la condici�n similar a la condici�n seco-mojado.
   call aguassomeras(modelo,navier,ap,tol,i,j,x,y,z,it,itt,nu,ma,sino,bcg,ten,newton,At,vn,vo,eval,ql,vib,vt,vb,est,imp,vit) 	      
   endif
  endif
  if ((maxval(eval).eq.0.0).or.(((it-app).gt.1).and.(tol.ne.1e-6))) then
  !O bien toda la soluci�n previa es superficial. 
  !No se calcula soluci�n subterr�nea y se cambia eval para que no indique posteriormente que toda el agua es subterr�nea en la siguiente
  !iteraci�n conjunta.
  !O bien se ha alcanzado convergencia para el modelo superficial (m�nimo de 2 iteraciones)
   do u=1,i
   eval(u)=sqrt(3.0_8)
   enddo
  else
  !La solucion previa es subterr�nea o superficial y subterr�nea
  !Siempre necesario usar la subrutina mallasubterranea previamente para obtener las alturas superficiales a aplicar en el contorno m�vil 
  !(que ir�n en vov(i)). Aqu� tambi�n se selecciona la malla para flujo subterr�neo, definida por el dominio superficial (la que malla que queda).   						 
   call mallasubterranea (i,j,z,zp,zzp,vt,eval,qx,qy,vov)
   if (minval(eval).eq.0.0) then
   !Hay soluci�n subterr�nea en al menos un elemento (si sol. subterr�nea s�lo en nodos aislados no se calcula con la ec. subterr�nea).
   call aguassubterranea(i,j,x,y,z,zp,zzp,nd,kix,kiy,ag,sino,bcg,At,vn,vov,eval,vb,ql,qb,qx,qy,vibv,vtv,velx,vely) 	   !*** A�adidos z,zp,vb
   endif
  endif
   if (((it-app).eq.1).and.(tol.ne.1e-6)) then
   !Se hace un m�nimo de 2 iteraciones conjuntas con Sw por paso de tiempo. Si no se hiciese esto se har� 1 iteraci�n si la diferencia es m�nima 
   !en la soluci�n superficial para el siguiente paso de tiempo lo cual es improbable (o tras la resoluci�n con 2D NS) y no se tendr�a la soluci�n
   !subterr�nea. Tampoco se tendr�a un buen procedimiento si se calculase. No se considerar�a la soluci�n subterr�nea obtenida en ese paso de 
   !tiempo en valores de velocidad.
   tol=1e-6
   write(6,*)'Se hace un minimo de 2 iteraciones conjuntas'
   endif	   
  enddo
  write(6,*)'Convergencia modelo conjunto en iteracion:', it
  if (sino.eq.'si')then
  write(6,*)'At=',Ata+At
  endif	 
endif 

!Arreglo de la soluci�n (interpolaci�n en nodos cuadr�ticos y c�lculo de variables conjuntas para el modelo conjunto).
!---------------------------------------------------------------------------------------------------------------------
if ((modelo.eq.'superficial').or.(modelo.eq.'conjunto')) then	  
!Para la pr�xima iteraci�n no-lineal ec. superficial en el siguiente paso de tiempo (transitorio) para el esquema impl�cito (ambos modelos) o para 
!la pr�xima iteraci�n no-lineal en el siguiente paso de tiempo para el esquema semi-impl�cito (modelo superficial), se considera la soluci�n 
!superficial obtenida (necesarios calados en ella):
do u=1,3*i
vit(u)=vt(u)
enddo
!Las alturas en los nodos esquina calculadas debajo difieren de las calculadas en 'vb' en que llevan las cotas del terreno en el subdominio 
!subterr�neo (modelo conjunto). As�, saldr� la cota de la l�mina de agua si hay calado y la cota del terreno si el calado es nulo.  
do u=1,i
vt(2*i+u)=vt(2*i+u)+z(u)
vx(u)=0.0_8
vy(u)=0.0_8
ht(u)=0.0_8
enddo
open(unit=1,file='C:\malla.txt',status='old')
!Operaciones hechas en los elementos del dominio superficial.
!Las alturas ir�n en la variable conjunta ht (se completar� m�s adelante en el caso del modelo conjunto).
!Las velocidades ir�n en las variables conjuntas vx,vy (se completar� m�s adelante en el caso del modelo conjunto). 
read(1,21)j
!Necesario leer j con los elementos superficiales (j ser� el n�mero de elementos subterr�neos si se ha usado el modelo conjunto).
read(1,'(A)')a
do u=1,j
 read(1,26)no(1),no(4),no(2),no(5),no(3),no(6) 
   do ui=1,3
   b=ui+1-sb(ui)
   c=ui+2-sb(ui+1)   
   !Interpolacion lineal postproceso de alturas y Manning en los nodos cuadr�ticos.
   vt(2*i+no(ui+3))=(vt(2*i+no(ui))+vt(2*i+no(b)))/2.0_8  
   ma(no(ui+3))=(ma(no(ui))+ma(no(b)))/2.0_8
   !Se escriben los valores de las variables para flujo superficial en las variables conjuntas.
   vx(no(ui))=vt(no(ui))	 
   vx(no(ui+3))=vt(no(ui+3))
   vy(no(ui))=vt(i+no(ui))
   vy(no(ui+3))=vt(i+no(ui+3))
   ht(no(ui))=vt(2*i+no(ui))
   ht(no(ui+3))=vt(2*i+no(ui+3))
   enddo
enddo
close(1)
endif
if ((modelo.eq.'subterraneo').or.(modelo.eq.'conjunto')) then
!Para la pr�xima iteraci�n no-lineal ec. subterr�nea en el siguiente paso de tiempo (transitorio) para el esquema impl�cito (ambos modelos), se 
!considera la soluci�n subterr�nea obtenida (necesarias alturas en ella):
do u=1,i
vitv(u)=vtv(u)
enddo
open(unit=3,file='C:\mallasub.txt',status='old')
!Operaciones hechas en los elementos del dominio subterr�neo.
read(3,21)j
read(3,'(A)')a
do u=1,j
read(3,26)no(1),no(4),no(2),no(5),no(3),no(6) 
 do ui=1,3
 b=ui+1-sb(ui)
 c=ui+2-sb(ui+1)
 !Interpolaci�n lineal postproceso de la velocidad subterr�nea, la altura (nivel fre�tico) y las propiedades en los nodos cuadr�ticos.
 vtv(no(ui+3))=(vtv(no(ui))+vtv(no(b)))/2.0_8
 kix(no(ui+3))=(kix(no(ui))+kix(no(b)))/2.0_8
 kiy(no(ui+3))=(kiy(no(ui))+kiy(no(b)))/2.0_8
 ag(no(ui+3))=(ag(no(ui))+ag(no(b)))/2.0_8
 nd(no(ui+3))=(nd(no(ui))+nd(no(b)))/2.0_8
 velx(no(ui+3))=(velx(no(ui))+velx(no(b)))/2.0_8
 vely(no(ui+3))=(vely(no(ui))+vely(no(b)))/2.0_8
 !Se completan las variables conjuntas escribiendo ahora en ellas los valores de las variables para flujo subterr�neo.
 vx(no(ui))=velx(no(ui))	 
 vx(no(ui+3))=velx(no(ui+3))
 vy(no(ui))=vely(no(ui))
 vy(no(ui+3))=vely(no(ui+3))
 if (zp(no(ui)).eq.zzp(no(ui)))then	                                !*** A�adida esta l�nea
 ht(no(ui))=vtv(no(ui))
 endif                                                              !*** A�adida esta l�nea
 if ((zp(no(ui)).eq.zzp(no(ui))).or.(zp(no(b)).eq.zzp(no(b)))) then !*** A�adida esta l�nea
 ht(no(ui+3))=vtv(no(ui+3))
 endif                                                              !*** A�adida esta l�nea
 
 

  !Se dar�n valores de calado en ciertos nodos cuadr�ticos que no tienen ecuaci�n superficial     !c�  
  !Para las condiciones puestas, no(ui) pertenecer� al dominio superficial, no(b) al subterr�neo y no(ui+3) ser� el n�do cuadr�tico entre ellos.   
  !Con valores de nivel fre�tico:
  !if ((vtv(no(ui+3)).gt.z(no(ui+3))).and.(vt(2*i+no(ui)).gt.z(no(ui))).and.(vt(2*i+no(b)).eq.z(no(b)))) then
  !vt(2*i+no(ui+3))=vtv(no(ui+3))
  !elseif ((vtv(no(ui+3)).gt.z(no(ui+3))).and.(vt(2*i+no(ui)).eq.z(no(ui))).and.(vt(2*i+no(b)).gt.z(no(b)))) then
  !vt(2*i+no(ui+3))=vtv(no(ui+3))
  !endif
  !Con valores de altura de nodos cercanos del dominio superficial (altura supuesta constante en el flujo superficial):
  if ((vt(2*i+no(ui)).gt.z(no(ui+3))).and.(vt(2*i+no(ui)).gt.z(no(ui))).and.(vt(2*i+no(b)).eq.z(no(b)))) then
  vt(2*i+no(ui+3))=vt(2*i+no(ui))
  !Para que la variable conjunta de la altura (ht) coincida con la variable altura superficial o cota del terreno (vt)
  !(con afecci�n en caso de aplicar el modelo conjunto): ht(no(ui+3))=vt(2*i+no(ui+3)) 
  elseif ((vt(2*i+no(b)).gt.z(no(ui+3))).and.(vt(2*i+no(ui)).eq.z(no(ui))).and.(vt(2*i+no(b)).gt.z(no(b)))) then
  vt(2*i+no(ui+3))=vt(2*i+no(b))
  !Para que la variable conjunta de la altura (ht) coincida con la variable altura superficial o cota del terreno (vt): ht(no(ui+3))=vt(2*i+no(ui+3))
  endif
 enddo
enddo
close(3)
endif

!Escritura de la soluci�n en fichero:
!------------------------------------
!Se escribe en el fichero solucionaguassomeras-subterraneo con formatos preparados para el programa Tecplot.	  
write(6,*)' '
if (sino.eq.'si') then
!Escritura para el esquema impl�cito con caso transitorio o para el esquema	semi-impl�cito (permite posterior animaci�n).
Ata=Ata+At
itt=itt+1
 !Se decide si se imprime la soluci�n (para algunos incrementos de tiempo podr�a no hacerse).
 if (ut.eq.frec*ni) then
 !Se pasan n�meros a formato caracter y se concatenan cadenas de caracteres. 
 write(Atee,'(i10)')idint(Ata)
 Atai='"solucion para t='//Atee
 Atai=Atai(1:27)//'s"'
 endif					 
else
!Escritura de una sola soluci�n para el esquema impl�cito con caso estacionario.
Atai='"solucion estacionaria unica"'
endif
if ((ut.eq.frec*ni).or.(sino.eq.'no')) then	 	  
write(8,22)Atai
endif					 

!Escritura de las variables representativas (con doble precisi�n) para cada modelo.
if (modelo.eq.'superficial') then
 !Modelo superficial. Para el siguiente paso de tiempo ec. superficial (transitorio) para el esquema impl�cito o para el siguiente paso de tiempo 
 !para el esquema semi-impl�cito, se considera como soluci�n en el instante anterior la obtenida (necesarias alturas en ella):
 do u=1,3*i
 vib(u)=vt(u) 
 enddo
 !Se decide si se imprime la soluci�n (para algunos incrementos de tiempo podr�a no hacerse).
 if ((ut.eq.frec*ni).or.(sino.eq.'no')) then	 
 open(unit=4,file='C:\mallainicial.txt',status='old')
 !Impresi�n de valores de coordenadas x e y en metros, de altura de la l�mina superficial-terreno en metros, de calado (incluido en NS-2D) en metros, 
 !de velocidad superficial (en direcciones x e y) en m/s, del m�dulo de la velocidad, de tensi�n (en direcciones x e y) en kilopascales, 
 !de tensi�n tangencial xy en kilopascales (kPa=1000Pa), y de vorticidad en 1/s.
 read(4,21)j
 write(8,23)i,4*j
 write(8,'(A)')'DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE )'
 do u=1,i
 te(u)=te(u)/1000.0_8
 te(i+u)=te(i+u)/1000.0_8
 te(2*i+u)=te(2*i+u)/1000.0_8
 write(8,38)x(u),y(u),vt(2*i+u),vt(2*i+u)-z(u),vt(u),vt(i+u),sqrt(vt(u)**2.0_8+vt(i+u)**2.0_8),te(u),te(i+u),te(2*i+u),vor(u)
 enddo
 endif					 
elseif (modelo.eq.'subterraneo') then
 !Modelo subterr�neo. Para el siguiente paso de tiempo ec. subterr�nea (transitorio) se considera como soluci�n en el instante anterior la obtenida
 !(necesarias alturas en ella):
 do u=1,i
 vibv(u)=vtv(u) 
 enddo
 !Se decide si se imprime la soluci�n (para algunos incrementos de tiempo podr�a no hacerse).
 if ((ut.eq.frec*ni).or.(sino.eq.'no')) then	 
 open(unit=4,file='C:\mallasubinicial.txt',status='old')
 !Impresi�n de valores de coordenadas x e y en metros, de nivel fre�tico en metros, de espesor fre�tico en metros, de velocidad de Darcy 
 !(en direcciones x e y) en m/s, y del m�dulo de la velocidad.
 read(4,21)j
 write(8,23)i,4*j
 write(8,'(A)')'DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE )'
 do u=1,i
 write(8,39)x(u),y(u),vtv(u),vtv(u)-zp(u),velx(u),vely(u),sqrt(velx(u)**2.0_8+vely(u)**2.0_8)
 enddo
 endif	                 
else
 !Modelo conjunto. Para el siguiente paso de tiempo ec. superficial y subterr�nea (transitorio) se considera como soluci�n en el instante anterior la 
 !soluci�n conjunta obtenida (variables conjuntas). Tambi�n es necesaria una soluci�n con las alturas en ella:
 do u=1,i
 vib(u)=vx(u)
 vib(i+u)=vy(u)
 !*** Primera opci�n
 if (zp(u).ne.zzp(u))then	!*** A�adida esta l�nea
 vib(2*i+u)=vt(2*i+u)		!*** A�adida esta l�nea
 vibv(u)=vtv(u)				!*** A�adida esta l�nea
 else						!*** A�adida esta l�nea
 vib(2*i+u)=ht(u)
 vibv(u)=ht(u)
 endif					    !*** A�adida esta l�nea
 !*** Otra opci�n (parece mejor dd no hay sol temporal de un modelo se da el terreno)
 !*** vib(2*i+u)=vt(2*i+u)
 !*** if (vtv(u).ne.0.0_8)then
 !*** vibv(u)=vtv(u)
 !*** else
 !*** vibv(u)=z(u)
 !*** endif
 !*** if (zp(u).ne.zzp(u))then
 !*** ht(u)=vt(2*i+u)
 !*** endif
 enddo
 if ((ut.eq.frec*ni).or.(sino.eq.'no')) then
 !Se decide si se imprime la soluci�n (para algunos incrementos de tiempo podr�a no hacerse).	 
 open(unit=4,file='C:\mallainicial.txt',status='old')
 !S�lo aqu� se utilizan de las variables conjuntas ht, vx y vy.
 !Impresi�n de valores de coordenadas x e y en metros, de altura de la l�mina superficial-nivel fre�tico (variable conjunta) en metros, de altura de la 
 !l�mina superficial-terreno en metros, de calado (incluido en NS-2D) en metros, de espesor fre�tico en metros, de velocidad superficial-velocidad de 
 !Darcy (en direcciones x e y, variables conjunta) en m/s, y del m�dulo de la velocidad.
 read(4,21)j	
 write(8,23)i,4*j
 write(8,'(A)')'DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE )'
 do u=1,i
 !Hay que tener en cuenta que el modelo conjunto no calcula flujo subterr�neo bajo el superifical.
 !Por defecto, siempre habr� calado nulo donde no se calcula flujo superficial (vt-z=vt original). Pero no ocurrir� lo mismo con el espesor fre�tico
 !donde no se calcula flujo subterr�neo. Aunque podr�a darse z(u)-zp(u), para ser consecuente se dar� espesor fre�tico nulo.
 !Ello permitir� reconocer hasta donde resuelve cada sub-modelo al representar estas variables.
 if (vtv(u).ne.0.0_8)then
 !En caso de calcularse flujo subterr�neo se imprimen todas las variables (el calado ser� nulo).
 write(8,40)x(u),y(u),ht(u),vt(2*i+u),vt(2*i+u)-z(u),vtv(u)-zp(u),vx(u),vy(u),sqrt(vx(u)**2.0_8+vy(u)**2.0_8)
 else
 !En caso de no calcularse flujo subterr�neo se da espesor fre�tico nulo. 
 write(8,40)x(u),y(u),ht(u),vt(2*i+u),vt(2*i+u)-z(u),0.0_8,vx(u),vy(u),sqrt(vx(u)**2.0_8+vy(u)**2.0_8)
 endif
 enddo		
 endif  	                 
endif

!Escritura de la malla.
if ((ut.eq.frec*ni).or.(sino.eq.'no')) then	 
 ni=ni+1					 
 write(8,'(A)')' '
 read(4,'(A)')a
 do u=1,j
 read(4,26)no(1),no(4),no(2),no(5),no(3),no(6)
 !Para tecplot cada elemento ha de ser lineal por lo que cada elemento supone cuatro a representar. En total 4*j elementos.
 !Por tanto se apreciar� una malla equivalente al representarla.
 write(8,41) no(1),no(4),no(6)
 write(8,41) no(4),no(2),no(5)
 write(8,41) no(5),no(3),no(6)
 write(8,41) no(4),no(5),no(6)
 enddo
 close(4)
 write(6,*)'Solucion escrita en fichero solucionaguassomeras-subterraneo'
endif                   
enddo

!Se decide si se termina la modelizaci�n:
!----------------------------------------
if (sino.eq.'si') then
write(6,*)'Quieres obtener la simulacion para un tiempo mayor? (si/no)'
read(5,*)sino
  if (sino.eq.'si') then
  nt=0
  ni=1				    
  endif
endif
enddo

!Escritura de propiedades (constantes para todos los incrementos de tiempo):
!---------------------------------------------------------------------------
!Se evita escribir datos innecesario en caso de utilizar el caso transitorio con el esquema impl�cito o el esquema semi-impl�cito. 
!Se escribe en el fichero solucionpropiedades con formatos preparados para el programa Tecplot.
if (modelo.eq.'superficial') then
 open(unit=4,file='C:\mallainicial.txt',status='old')
 read(4,21)j
 !Valores de coordenadas en metros (con cota del terreno), y de Manning,
 write(10,'(A)')'TITLE = "Propiedades para aguassomeras"'
 write(10,'(A)')'VARIABLES =  X, Y, Z, n'
 write(10,24)i,4*j
 write(10,'(A)')'DT=(DOUBLE DOUBLE DOUBLE DOUBLE )'
 do u=1,i
 write(10,48)x(u),y(u),z(u),ma(u)
 enddo 
elseif (modelo.eq.'subterraneo') then
 open(unit=4,file='C:\mallasubinicial.txt',status='old')
 read(4,21)j
 write(10,'(A)')'TITLE = "Propiedades para flujo subterraneo"'
 !Valores de coordenadas en metros (con cota del sustrato), de conductividad (en direcciones x e y) en m/s, de �ngulo de anisotrop�a en grados, 
 !y de porosidad (adimensional).
 write(10,'(A)')'VARIABLES =  X, Y, ZP, Kix, Kiy, Ang�, nd'
 write(10,24)i,4*j
 write(10,'(A)')'DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE )'
 do u=1,i
 write(10,49)x(u),y(u),zp(u),kix(u),kiy(u),ag(u)*180.0_8/3.14159_8,nd(u)
 enddo
else
 open(unit=4,file='C:\mallainicial.txt',status='old')
 read(4,21)j
 write(10,'(A)')'TITLE = "Propiedades para aguassomeras y subterraneo"'
 !Valores de coordenadas en metros (con cota del terreno y cota del sustrato), de conductividad (en direcciones x e y) en m/s, de �ngulo 
 !de anisotrop�a en grados sexagesimales, de porosidad (adimensional) y de Manning.
 write(10,'(A)')'VARIABLES =  X, Y, Z, ZP, Kix, Kiy, Ang�, nd, n'
 write(10,24)i,4*j
 write(10,'(A)')'DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE )'
 do u=1,i
 write(10,50)x(u),y(u),z(u),zp(u),kix(u),kiy(u),ag(u)*180.0_8/3.14159_8,nd(u),ma(u)
 enddo
endif
!Escritura de la malla
write(10,'(A)')' '
read(4,'(A)')a
do u=1,j
read(4,26)no(1),no(4),no(2),no(5),no(3),no(6)
write(10,41) no(1),no(4),no(6) 
write(10,41) no(4),no(2),no(5) 
write(10,41) no(5),no(3),no(6) 
write(10,41) no(4),no(5),no(6) 
enddo
close(4)

!Eliminaci�n de los ficheros generados:
!--------------------------------------
if ((modelo.eq.'superficial').or.(modelo.eq.'conjunto')) then
open(unit=1,file='C:\malla.txt',status='old')
close(1,status = 'delete')
endif
if ((modelo.eq.'subterraneo').or.(modelo.eq.'conjunto')) then
open(unit=3,file='C:\mallasub.txt',status='old')
close(3,status = 'delete')
endif

deallocate(x,y,z,zp,zzp,vt,vit,vib,vb,vo,eval,vov,vtv,vitv,vibv,kix,kiy,ag,nd,ma,qx,qy,qb,ql,velx,vely,vx,vy,ht,evol,te,vor)
end

!----------------------------------------------------------------------------------------------------------------------------------------------------------
!Subrutinas de gesti�n para el c�lculo de la soluci�n de las ecuaciones.
!----------------------------------------------------------------------------------------------------------------------------------------------------------
!Subrutina AGUASSOMERAS (tres ecuaciones, 2 din�micas y 1 de continuidad). En esta subrutina se calculan las velocidades, la altura de la l�mina de agua. 
!Se construye el sistema con la matriz de rigidez (de cajas en las que ir�n A, B,...). Se calculan las alturas (respecto a z=0) en nodos esquina y las 
!velocidades en todos los nodos. Se calculan los calados postproceso. 
!Se pueden resolver por el m�todo de Picard la ecuaci�n de aguas someras o la ecuaci�n de N-S 2D.
!Se puede utilizar el m�todo de Newton para resolver la ecuaci�n de aguas someras, en cuyo caso se aplicar�n 5 iteraciones previas de Picard.
!S�lo se har� una iteraci�n no-lineal y se construir�n todas las cajas cada vez (la posible afecci�n de la condici�n seco-mojado obliga 
!a ello aunque las cajas no contengan coeficientes no lineales).  
!El vector inicial (it=0, itt=0) para resolver la no linealidad estar� en vt y vb (aprox inicial en estacionario, CI u otro valor en transitorio) y tendr�
!los valores Dirichlet prescritos (variables del sistema) de la altura de la l�mina de agua y de velocidad para el c�lculo de las integrales de contorno.
!Al realizarse la primera iteraci�n de Picard (it=0, itt=0) se pueden resolver algunas iteraciones de las ecuaciones de N-S 2D al resolver 
!la ecuaci�n de aguas someras. As� se calcula un buen vector inicial para resolver en la siguiente iteraci�n las ecuaciones de aguas someras.
!Siempre habr� un calado, altura y velocidad nulos donde no se calcula flujo superficial (en vt y vb).
!'tol' gestiona la parada del bucle cuando dos soluciones consecutivas son suficientemente similares,
!'vn,ve,ven,vo,vu' son contadores de nodos cuadr�ticos, lineales o con condici�n de contorno,
!'veac,veoc,veuc' son vectores que contienen las integrales de contorno,'vef' es el t�rmino fuente de caudal por lluvia para flujo superficial,
!'fx,fy' son los vectores donde se calculan las integrales de la pendiente de fricci�n (cuantifican la tensi�n de fondo), 
!'fz' es un vector del mismo tipo que se genera al aplicar el m�todo de estabilizaci�n,
!'sa' es un vector donde se almacenan los coeficientes no nulos de la matriz del sistema (formato MSR), 
!'ita' es el vector puntero con la posici�n de estos coeficientes (formato MSR),
!'cia,ca' son inicialmente an�logos a 'ita,sa' al ser usados para guardar por separado las matrices de masa. Si se resuelve con 
!precondicionador diagonal (subrutina gradientesbiconjugados para formato MSR) 'ita,sa' llevar�n la matriz del sistema a resolver y 'cia,ca' son desechados.
!En otro caso se dimensiona 'cja', se copiar�n 'ita,sa' a 'ca,cia,cja' (subrutina dslubc para formato CSC) y 'ita,sa son desechados. As�,
!'ca' ser� el vector donde se almacenan los coeficientes no nulos de la matriz del sistema (formato CSC), 
!'cia,cja' ser�n los vectores punteros con la posici�n de estos coeficientes (formato CSC),
!'ndim' es el n� de coeficientes no nulos que hay en la matriz del sistema almacenados en los vectores 'ita,sa' (varia al formar el sistema a resolver),
!'inc' es la dimensi�n de la matriz del sistema cuadrada si se usasen elementos cuadr�ticos para alturas/calados y velocidades (n�ecuaciones*i), 
!'c' es la dimensi�n de la matriz del sistema cuadrada para resolver el sistema (en cada iteraci�n) tras reducir el orden con las CC y considerar los 
!elementos lineales para alturas/calados, 
!'ru,rb' son los vectores del residual y de la soluci�n para el m�todo de Newton,
!'vec,voc' forman el t�rmino independiente del sistema, 'vv,vt,vb' son los vectores donde se guarda la soluci�n del sistema (en cada iteraci�n), 
!'yn' permite el c�lculo de una soluci�n de aguas someras con calados similares sea como sea la cota del terreno,
!'nonzero' permite utilizar la soluci�n anterior para aplicar el m�todo de los gradientes biconjungados,
!'vth' lleva la soluci�n del calado para la aplicaci�n de la condici�n seco-mojado (subrutina nuevamalla).                                                                               
!-----------------------------------------------------------------------------------------------------------------------------------------------------------
subroutine aguassomeras(modelo,navier,ap,tol,i,j,x,y,z,it,itt,nu,ma,sino,bcg,ten,newton,At,vnin,vo,eval,ql,vib,vt,vb,est,imp,vit)
use allocatacion
integer*4, dimension(:),allocatable::vn
integer*4 i,j,u,uu,c,vnin(i),it,itt,inc,ndim,k,conv			 
real*8, dimension(:),allocatable::vbdin,vecdin,vv,vec,voc,vu,vef,fx,fy,fz,vth,veac,veoc,veuc,ru,rb		  	 
real*8 x(i),y(i),z(i),ma(i),vt(3*i),vb(3*i),vit(3*i),vib(3*i),tol,nu,vo(3*i),eval(i),At,ql(i),del									  		
logical nonzero 
character modelo*12,sino*2,newton*2,bcg*2,ten*2,navier*2,yn*2,est*2,imp*2,ap*2
!Aqu� s� hace un dimensionamiento din�mico propiamente dicho con 'vbdin,vecdin' (ser�a obligatorio definirlas as� si fuesen variables globales). 
!Otro ejemplo son las variables globales ita,sa (usadas para resolver con precondicionador diagonal) � cia,ca,cja (usadas para resolver con 
!precondicionador LU) cuya memoria se destruir� (seg�n el caso) en tiempo de ejecuci�n de esta subrutina (tras resolver con ellos).
allocate(vn(i),vv(3*i),vec(3*i),voc(3*i),vu(3*i),vef(i),fx(i),fy(i),fz(i),vth(i),veac(i),veoc(i),veuc(i),ru(3*i),rb(3*i))

30  format(2(3x,A1,'(',I5,')=',E15.8E2))	  

!Condici�n 'falsa' de calado normal. !c�
!---------------------------------------
!Si yn='si' se genera una soluci�n con calado aproximadamente id�ntico en todos los nodos independiente de la cota del terreno.
!S�lo tenida en cuenta para newton='no' y est='no' y con sentido si navier='no' (aguas someras).
!Si yn='no' se resuelven las ecuaciones apropiadas.
yn='no'  
!Para la aplicaci�n del esquema impl�cito o semi-impl�cito (el semi-impl�cito ser� de segundo orden y conllevar� a menor error).
if (imp.eq.'no') then
del=2.0_8
else
del=1.0_8
endif

!Selecci�n de en qu� iteraciones se utilizar� la soluci�n anterior para el m�todo de los gradientes biconjungados
!(aunque en caso de aplicar el m�todo el Newton esto no ser� efectivo)
!----------------------------------------------------------------------------------------------------------------
if (it.eq.0)then
nonzero=.false.	
else
nonzero=.true.
endif

!Inizializaci�n de variables:
!----------------------------
inc=3*i
do u=1,3*i
vv(u)=vt(u)
ru(u)=0.0_8
rb(u)=vb(u)
voc(u)=0.0_8
vu(u)=vo(u)
enddo
do u=1,i
fx(u)=0.0_8
fy(u)=0.0_8
fz(u)=0.0_8
vef(u)=0.0_8
veac(u)=0.0_8
veoc(u)=0.0_8
veuc(u)=0.0_8
enddo
!El vn original es ahora vnin. El nuevo vn ser� diferente si el dominio es menor al de toda la malla.
!En la ecuaci�n de continuidad se eliminar�n las ecuaciones relativas a todos los nodos cuadr�ticos y no s�lo los del dominio superficial.
!Con vn se eliminar�n los de este dominio y el resto se eliminar�n con las CC al no considerar el resto de la malla.
do u=1,i
 if ((vnin(u).eq.1).and.(eval(u).eq.0.0)) then
 !Nodos cuadr�ticos y dominio superficial
 vn(u)=1
 else
 vn(u)=0
 endif
enddo
!Se dimensionan ita y sa (tambi�n cia y ca), de acuerdo al n�mero total de coeficientes que se generar�n con las matrices elementales.
!Adem�s ita ya lleva referenciados el n�mero de coeficientes que habr� por fila (en sus primeras 3*i+1 componentes) antes de dar el formato MSR.
!Cuando se calcule una matriz se almacenar�n sus coeficientes directamente en estos vectores.
call dimvectas (i,j,newton,sino,vnin,est)
ndim=ita(inc+1)-1

!En caso de resolver de forma transitoria con incrementos de tiempo (esquema impl�cito considerando las matrices de masa o esquema semi-impl�cito).
!--------------------------------------------------------------------------------------------------------------------------------------------------
!Si se aplica el modelo conjunto y tiene afecci�n la condici�n similar a la condici�n seco-mojado, se debe tener en cuenta que las condiciones
!de contorno que se impondr�n en el contorno m�vil (que no coincidir�n con los valores de la soluci�n temporal previa y que son 
!generadas por el modelo subterr�neo) se corresponden con un tiempo igual al tiempo para el que se busca la soluci�n.
if (sino.eq.'si') then
call timeasu (i,j,x,y,At)
!Para formulaci�n BG las matrices M (timeasu) valen para las ecuaciones de aguas someras y de N-S 2D. 
!Para formulaci�n estabilizada esto no es as� para las nuevas matrices de masa estabilizadas porque se usan par�metros distintos. 
!De todos modos no funciona bien la consideraci�n de matrices de masa estabilizadas y de momento no se usan (l�neas comentadas debajo).
 if (navier.eq.'si') then 
  !if (est.eq.'si')then
  !call timeasupgns (i,j,x,y,vv,At)
  !else
  continue
  !endif
 else
  call timeasd (i,j,x,y,At)
  !if (est.eq.'si')then
  !call timeasupgas (i,j,x,y,vv,At)
  !endif
 endif
!Se guarda la matriz de rigidez con las matrices de masa.
do u=1,ndim
cia(u)=ita(u)
ca(u)=sa(u)
enddo
!Se ordenan los coeficientes (alguno podr�a ser nulo), se ensamblan y se eliminan los coeficientes nulos. Resumiendo, se da el formato MSR.
!S�lo as� es posible calcular el producto de la matriz de rigidez (con las matrices de masa) por la soluci�n en el instante anterior.
call orden (inc,k)
do u=1,3*i
 do uu=ita(u),ita(u+1)-1
 voc(u)=voc(u)+sa(uu)*vib(ita(uu))
 enddo
 voc(u)=voc(u)+sa(u)*vib(u)
enddo
!Se vuelve a tener la matriz de rigidez en 'ita,sa' sin ordenar. As�, ser� posible sumar directamente los coeficientes correspondientes a otras
!matrices que est�n en la misma caja de la matriz (por ejemplo A y M) sin buscar donde deben ser almacernardos (dado que el orden con que se toman 
!los elementos de la malla ser� el mismo). 
do u=1,ndim
ita(u)=cia(u)
sa(u)=ca(u)
enddo
endif 	   

!Se calculan las matrices A y B de las ecuaciones din�micas. Valen para las ecuaciones de aguas someras y de N-S 2D:
!-------------------------------------------------------------------------------------------------------------------
call cajasab (i,j,x,y,nu,del)  

!En la primera iteraci�n se sustituyen las CC en el vector inicial para calcular las integrales de contorno.
!-----------------------------------------------------------------------------------------------------------
!En el resto de iteraciones se utilizar� la soluci�n y los valores ser�n los mismos a los de la CC al haber reducido el tama�o del sistema con ellas.
if ((it.eq.0).and.(itt.eq.0).and.(navier.eq.'si'))then
!En caso de (navier.eq.'no') esto ya se ha hecho en la subrutina mallaaguassomeras. Se recuerda que en este caso se hace la selecci�n inicial 
!y es necesaria esa subrutina. Adem�s, en este caso no ser�a necesario hacer esto si la CI o aproximaci�n inicial tiene los valores de las CC.
 do u=1,3*i
  if (vo(u).ne.sqrt(2.0_8)) then 
   if (u.gt.2*i)then
   vv(u)=vo(u)-z(u-2*i)
   vb(u)=vo(u)
   else
   vv(u)=vo(u)
   endif  
  endif
 enddo
endif
 
!Otras matrices y vectores para las ecuaciones de Navier-Stokes 2D:
!------------------------------------------------------------------
if (navier.eq.'si') then        
 !C�lculo de las integrales de contorno	para el caso impl�cito.
 call vectorcontornopresiones(i,j,x,y,vv,vb,veac,veoc,veuc,nu,ten) 
 if (imp.eq.'no')then
 !Suma del c�lculo de otras integrales de contorno para el caso semi-impl�cito.
 call vectorcontornopresiones(i,j,x,y,vit,vib,veac,veoc,veuc,nu,ten)
 endif
  !C�lculo de las matrices Bt de la ecuaci�n de continuidad
  call cajasbt (i,j,x,y,del)
  !Suma de las integrales de contorno al t�rmino independiente. 
  do u=1,i
  vec(u)=del*voc(u)+veac(u)/del
  vec(i+u)=del*voc(i+u)+veoc(u)/del
  vec(2*i+u)=voc(2*i+u)
  enddo
   !Suma de una matriz (m�s completa) por un vector al t�rmino independiente para el caso semi-impl�cito.
   if (imp.eq.'no')then
    do u=1,ndim
    cia(u)=ita(u)
    ca(u)=sa(u)
    enddo
    call matriznolineal(i,j,x,y,vit,del)
    if (est.eq.'si')then
    call cajasupgns(i,j,x,y,vit,nu,del)
    endif
    call orden (inc,k)
    do u=1,3*i
     do uu=ita(u),ita(u+1)-1
     vec(u)=vec(u)-sa(uu)*vib(ita(uu))
     enddo
     vec(u)=vec(u)-sa(u)*vib(u)
    enddo
    do u=1,ndim
    ita(u)=cia(u)
    sa(u)=ca(u)
    enddo
   endif
  !C�lculo de las matrices estabilizadas.
  if (est.eq.'si')then
  call cajasupgns(i,j,x,y,vv,nu,del)
  endif    

!Otras matrices y vectores para las ecuaciones de aguas someras:
!---------------------------------------------------------------
else
 !C�lculo de las integrales de contorno	y las integrales de la pendiente de fricci�n para el caso impl�cito.
 call vectorcontornopresiones(i,j,x,y,vv,vb,veac,veoc,veuc,nu,ten)
 call f(i,j,x,y,z,ma,yn,vv,fx,fy)
  !Incluyen las integrales de la pendiente de fricci�n estabilizadas.
  if (est.eq.'si')then
  call fsupg(i,j,x,y,ma,vv,fx,fy,fz)
  endif
 if (imp.eq.'no')then
 !Suma del c�lculo de otras integrales de contorno y otras integrales de la pendiente de fricci�n para el caso semi-impl�cito.
 call vectorcontornopresiones(i,j,x,y,vit,vib,veac,veoc,veuc,nu,ten)
 call f(i,j,x,y,z,ma,yn,vit,fx,fy)
  if (est.eq.'si')then
  call fsupg(i,j,x,y,ma,vit,fx,fy,fz)
  endif
 endif
 !C�lculo de las integrales correspondientes a la lluvia.
 call vectorcontornofuente(i,j,x,y,vef,ql)	
  !Suma de las integrales �nteriores al t�rmino independiente.
  do u=1,i
  vec(u)=del*voc(u)-fx(u)*9.81_8/del+veac(u)/del
  vec(i+u)=del*voc(i+u)-fy(u)*9.81_8/del+veoc(u)/del
   if (est.eq.'si') then
   vec(2*i+u)=del*voc(2*i+u)-fz(u)*9.81_8/del-veuc(u)/del
   else
   vec(2*i+u)=del*voc(2*i+u)-fz(u)*9.81_8/del-veuc(u)/del+vef(u)
   endif
  enddo  
   !Suma de una matriz (m�s completa) por un vector al t�rmino independiente para el caso semi-impl�cito.
   if (imp.eq.'no')then
    do u=1,ndim
    cia(u)=ita(u)
    ca(u)=sa(u)
    enddo
    call cajasde (i,j,x,y,vit,del)
    call matriznolineal(i,j,x,y,vit,del)
    if (est.eq.'si')then
    call cajasupgas(i,j,x,y,vit,nu,del)
    endif
    call orden (inc,k)
    do u=1,3*i
     do uu=ita(u),ita(u+1)-1
     vec(u)=vec(u)-sa(uu)*vib(ita(uu))
     enddo
     vec(u)=vec(u)-sa(u)*vib(u)
    enddo
    do u=1,ndim
    ita(u)=cia(u)
    sa(u)=ca(u)
    enddo
   endif  
  !C�lculo de las matrices D y E de la ecuaci�n de continuidad.
  call cajasde (i,j,x,y,vv,del)
  !C�lculo de las matrices estabilizadas.
  if (est.eq.'si')then
  call cajasupgas(i,j,x,y,vv,nu,del)
  endif  
endif

!Se calculan las matrices C de las ecuaciones din�micas. Valen para las ecuaciones de aguas someras y de N-S 2D.
!---------------------------------------------------------------------------------------------------------------
call matriznolineal(i,j,x,y,vv,del)	  

!Se ha aplicado un esquema semi-impl�cito determinado (casi de Crank-Nicolson)
!Con vit se tiene otra versi�n seg�n el libro de Quarteroni y Valli (chap. 13). De programarla, s�lo llamar�a una vez a cada caja con vit 
!y me traer�a hasta aqu� el trozo "donde se suma una matriz (mas completa) por un vector al t�rmino independiente". 
!En caso de hacerlo as�, no podr�a tener programado a la vez el impl�cito.

!M�todo de Picard:
!-----------------
if (((it.lt.5).and.(newton.eq.'si')).or.(newton.eq.'no')) then
 !Se imponen las CC sobre el sistema y se forma la matriz de dimensi�n (2*i + n�nodos esquina)-n� de CC para en el dominio;
 !Se reduce tambi�n en los nodos donde eval=sqrt(3)	si la condici�n seco-mojado (o la condici�n similar a la condici�n seco-mojado) tiene afecci�n.
 call reducciondelsistema(i,c,vec,vu,vn,vb,inc,bcg,ndim)	  												
 !Resolucion mediante gradientes biconjugados. Se utilizar� la soluci�n anterior (nonzero=.true.) para acelerar el m�todo.																		
 allocate (vbdin(c),vecdin(c))	 
 do u=1,c			
 vecdin(u)=vec(u)
  if (nonzero)then
  vbdin(u)=vb(u)  
  else
  vbdin(u)=0.0_8
  endif
 enddo
  !Uso de las subrutinas (se�aladas en Press et al.) que calculan con un precondicionador diagonal.
  conv=1
  u=0  
  if (bcg.eq.'no') then
    do while (conv.eq.1)
    call gradientesbiconjugados(vecdin,c,vbdin,nonzero,conv)
	nonzero=.true.
    u=u+1
    if ((u.eq.3).and.(conv.eq.1)) then
	!Se calcula yendo de C en C iteraciones con PCGB y haci�ndolo hasta 3 veces.
	!Se para el programa al llegar a tres veces y se le pide que est� atento a si hay problemas 
	!a partir de ah� y que escoja otro m�todo si es as�.
     if (imp.eq.'si') then
	 write(6,*) 'No es posible obtener solucion del sistema en iteracion no lineal',it
	 else
	 write(6,*) 'No es posible obtener solucion del sistema'
	 endif
	 write(6,*) 'Pulsa ENTER. Si hay problemas a partir de aqui escoge otro metodo'
	 read(5,*)
    exit
    endif
    enddo 
  !Uso de las subrutinas (incluidas en la librer�a SLATEC) que calculan con un precondicionador LU aproximado.
  !All� se definen 'nu' como el n� de coeficientes no nulos en la parte triangular superior de la matriz de rigidez, 'nl' el n� en la parte inferior 
  !y 'nelt' el n�mero total de coeficientes no nulos. Se usan las dimensiones nl+nu+8*n para el vector de coeficientes enteros y nl+nu+4*n+2 para el 
  !vector de coeficientes reales. Se cumplir� que nu+nl = nelt-n = ndim-1-n = ndim-1-c. 
  else                                  	
    do while (conv.eq.1)                                  	
    call dslubc (vecdin,c,vbdin,ndim+7*c,ndim+3*c+12,conv)	
    u=u+1
    if ((u.eq.3).and.(conv.eq.1)) then
	!Se calcula yendo de C en C iteraciones con PCGB y haci�ndolo hasta 3 veces.
     if (imp.eq.'si') then
	 write(6,*) 'No es posible obtener solucion del sistema en iteracion no lineal',it
	 else
	 write(6,*) 'No es posible obtener solucion del sistema'
	 endif
	 write(6,*) 'Pulsa ENTER. Si hay problemas a partir de aqui escoge otro metodo'
	 read(5,*)
    exit
    endif
    enddo	  
  endif
  do u=1,c			
  vb(u)=vbdin(u)
  enddo
 !Arreglo del vector soluci�n introduciendo los valores prescritos (que incluye valores nulos de velocidad y alturas para la parte de la malla seca 
 !en nodos cuadr�ticos y lineales si tiene afecci�n la condici�n seco-mojado o la condici�n similar a la condici�n seco-mojado) 
 !y valores nulos de altura en los nodos cuadr�ticos (s�lo en los del trozo mojado de la malla en caso de afecci�n de la condici�n).	
 do u=1,3*i
  if (vu(u).ne.sqrt(2.0_8)) then 
  !Si la condici�n (vu.ne.) con mayor u es para u<<3*i se entra con c=3*i-1. Se sale con c=3*i y al bucle 'u' a�n le quedar�n iteraciones (no se entrar�).
  !si la condici�n (vu.ne.) con mayor u es para u=3*i-1 se entra con c=3*i-1. Se sale con c=3*i y queda una �ltima iteraci�n de 'u' (no se entrar�).
  !Si la condici�n (vu.ne.) con mayor u es para u=3*i	se entra con c=3*i-1. As�, no se opera el siguiente if, se saldr� con c=3*i.
  !Al ir por orden 'u' s�lo puede ser una unidad mayor a 'c' y 'vb' ya estar� ordenado.    
   if (u.le.c) then
   do uu=c,u,-1
   vb(uu+1)=vb(uu)			  
   enddo
   endif
   c=c+1
   vb(u)=vu(u)
  endif 
 enddo

!M�todo de Newton: s�lo se puede usar para formulaci�n BG, para ec. aguas someras, para esquema impl�cito y para it>4
!--------------------------------------------------------------------------------------------------------------------
else
 !Aplicaci�n del m�todo a la ec. de aguas someras del mismo modo que se aplica a la ec. de N-S 2D en Reddy y Gartling. 
 !Suma de la matriz de rigidez por la soluci�n anterior del problema (rb=vb) y del antiguo t�rmino independiente cambiado de signo 
 !al nuevo t�rmino independiente.
 do u=1,ndim
 cia(u)=ita(u)
 ca(u)=sa(u)
 enddo
 call orden (inc,k)
 do u=1,3*i
  do uu=ita(u),ita(u+1)-1
  ru(u)=ru(u)+sa(uu)*rb(ita(uu))	   
  enddo
  ru(u)=ru(u)+sa(u)*rb(u)-vec(u)
 enddo
 do u=1,ndim	
 ita(u)=cia(u)
 sa(u)=ca(u)
 enddo 
 !Se calcula la matriz de las derivadas (matriz jacobiana), la nueva matriz del sistema.
 !Dentro de la subrutina se siguen los mismos pasos que en la subrutina reducciondelsistema, s�lo que en este caso se 
 !imponen CC nulas all� donde se conozca una CC.
 call jacob (i,j,x,y,vv,nu,ma,ru,vn,vu,vb,c,ten,bcg,ndim)
 allocate (vbdin(c),vecdin(c))	 
 !Resolucion mediante gradientes biconjugados. 
 !Mejor comportamiento inicializando a cero (nonzero=.false.) que guardar la �ltima soluci�n (la soluci�n se hara cero en la convergencia). 
 !Adem�s se evita almacenar una nueva variable ya que en vb no est� la �ltima soluci�n del sistema, (ha sido 
 !sobreescrito en la iteraci�n previa con vb=rb-vb y tiene la soluci�n anterior del problema).
 nonzero=.false.
 do u=1,c			
 vecdin(u)=ru(u)  
 vbdin(u)=0.0_8 
 enddo
  conv=1
  u=0
  if (bcg.eq.'no') then
    do while (conv.eq.1)
    call gradientesbiconjugados(vecdin,c,vbdin,nonzero,conv)
	nonzero=.true.
    u=u+1
    if ((u.eq.3).and.(conv.eq.1)) then
     if (imp.eq.'si') then
	 write(6,*) 'No es posible obtener solucion del sistema en iteracion no lineal',it
     else
	 write(6,*) 'No es posible obtener solucion del sistema'
	 endif
	 write(6,*) 'Pulsa ENTER. Si hay problemas a partir de aqui escoge otro metodo'
	 read(5,*)
	exit
    endif
    enddo
  else 
    do while (conv.eq.1)
	call dslubc (vecdin,c,vbdin,ndim+7*c,ndim+3*c+12,conv)
	u=u+1
    if ((u.eq.3).and.(conv.eq.1)) then
     if (imp.eq.'si') then
	 write(6,*) 'No es posible obtener solucion del sistema en iteracion no lineal',it
     else
	 write(6,*) 'No es posible obtener solucion del sistema'
	 endif
	 write(6,*) 'Pulsa ENTER. Si hay problemas a partir de aqui escoge otro metodo'
	 read(5,*)
	exit
    endif
    enddo	   
  endif
  !La soluci�n del sistema tendr� todos los valores m�s pr�ximos a cero tanto mayor n�mero de iteraciones se hagan.
  do u=1,c			
  vb(u)=vbdin(u)
  enddo
 !Arreglo del vector soluci�n del mismo modo que al terminar cada iteraci�n de Picard (los valores prescritos ahora tambi�n ser�n valores nulos).
 do u=1,3*i
  if (vu(u).ne.sqrt(2.0_8)) then
   if (u.le.c) then
   do uu=c,u,-1
   vb(uu+1)=vb(uu)
   enddo
   endif
   c=c+1
   vb(u)=0.0_8
  endif 
 enddo
 !Nueva soluci�n = soluci�n anterior - soluci�n del sistema
 do u=1,3*i
 vb(u)=rb(u)-vb(u)
 enddo
endif

!Postproceso tras la resoluci�n del sistema:
!-------------------------------------------
!Se calcula el calado en los nodos cuadr�ticos.
!Da igual que los calados de la soluci�n de aguas someras sean menor que cero porque 
!con la condici�n seco-mojado se crear� una nueva malla que no considerar� estos valores.
do u=1,3*i
vt(u)=vb(u)
enddo
do u=1,i	  
if ((vnin(u).ne.1).and.(eval(u).eq.0.0)) then
vt(2*i+u)=vt(2*i+u)-z(u)		                         		   
endif
enddo
it=it+1

!Se muestra la soluci�n de velocidades y calados por pantalla (s�lo en los nodos del dominio superficial si la condici�n
!seco-mojado o la condici�n similar a la condici�n seco-mojado tiene afecci�n).	Impresi�n de resultados fraccionada, de 200 en 200
!l�neas (necesario definir 'ua' como entero).
!write(6,*)' '
!ua=1
!uu=0
!do u=1,i
! if (eval(u).eq.0.0) then								 
! write(6,30)'u',u,vt(u),'v',u,vt(u+i)	
! uu=uu+1
!  if (uu.eq.200*ua) then
!  write(6,*)'Mostrados arriba 200 resultados. Pulsa enter.'
!  read(5,*)
!  ua=ua+1
!  endif  
! endif
!enddo
!write(6,*)' '
!ua=1
!uu=0
!do u=1,i
! if ((vnin(u).ne.1).and.(eval(u).eq.0.0)) then		             
! write(6,30)'h',u,vt(u+2*i)
! uu=uu+1
!  if (uu.eq.200*ua) then
!  write(6,*)'Mostrados arriba 200 resultados. Pulsa enter.'
!  read(5,*)
!  ua=ua+1
!  endif
! endif
!enddo
	 
!Se muestra por pantalla cu�ndo se da el caso particular de tener la soluci�n del flujo de Stokes 2D, en el que no hay convecci�n (se da ya que
!la matriz convectiva no lineal es nula al introducir velocidades nulas y no es necesario resolver un problema no lineal). 
!if ((it.eq.1).and.(itt.eq.0).and.(navier.eq.'si')) then 
!write(6,*)' '
!write(6,*)'Primera iteracion del metodo de Picard con ecuaciones de NS 2D partiendo'
!write(6,*)'de un valor nulo. Equivale a la solucion del flujo de Stokes 2D.'	  
!endif
	  		  
!Control para la parada del bucle para la iteraci�n no lineal.
if (imp.eq.'si') then
write(6,*)'Error entre iteraciones no lineales: ',maxval(abs(vv-vt))
endif
!Cuando la diferencia entre cualquier componente del vector soluci�n y la obtenida en la iteraci�n anterior sea menor que la variable tol
!definida en la subrutina aguassomerassubt se detendr�n las iteraciones (se sobreescribe aqu� tol y as� no se entrar� en esta subrutina).
if (maxval(abs(vv-vt)).lt.tol) then
tol=1.0
endif

if (bcg.eq.'no') then
deallocate(vbdin,vecdin,ita,sa)
else
deallocate(vbdin,vecdin,cia,ca,cja) 
endif

!Condici�n seco-mojado � condici�n similar a la condici�n seco-mojado:
!---------------------------------------------------------------------
!La condici�n seco-mojado deja la soluci�n preparada para otra iteraci�n del modelo superficial (da valores, selecciona dominio y da CC en contorno m�vil)
!Es para el modelo superficial.
!La condici�n similar deja la soluci�n preparada para el modelo subterr�neo (da valores). Antes de resolver el subterr�neo se usa una subrutina previa
!(que selecciona dominio y da CC). Tras resolver el subterr�neo tambi�n ser� necesaria otra subrutina previa (que seleccione dominio y de CC). 
!Es para el modelo conjunto.
do u=1,i
vth(u)=vt(2*i+u)
enddo
!Si se busca la soluci�n de las ecuaciones de N-S 2D no se aplica esta condici�n. Se aplica si se busca la soluci�n de las ec. de aguas someras
!con o sin iteraciones previas de N-S 2D (siempre sucede al aplicar el modelo conjunto). 
!Por tanto siempre que se entra aqu�, se habr� evaluado la condici�n inicial o la aproximaci�n inicial en la subrutina aguassomerassubt 
if ((navier.eq.'no').or.((navier.eq.'si').and.(ap.eq.'si'))) then 
if (maxval(vth).le.0.0) then 
 !Para it=1 (ya es 1 en la tras la 1� iteraci�n) si se resuelve con toda la malla y salen todos los calados negativos se tendr� maxval(vth).lt.0.0.
 !Para it>1 se tendr� un dominio menor (con vth=0 fuera), y si en �l salen todos los calados negativos se tendr� maxval(vth).le.0.0.   
 write(6,*)'El flujo superficial ha desaparecido'
 do u=1,i
 eval(u)=sqrt(3.0_8)
 enddo
 if (modelo.eq.'superficial') then
  do u=1,3*i
  vt(u)=0.0_8
  enddo
  tol=1.0
 endif
elseif ((minval(vth).lt.0.0).or.(maxval(eval).ne.0.0))then
 !Si se ha considerado toda la malla y los valores de calado son positivos no se entra aqu�. 
 !El utilizar una condici�n que eval�e el dominio considerado (eval) adem�s de otra que eval�e el signo del calado (vth), permite entrar y a�adir 
 !elementos al nuevo dominio a�n cuando no se obtengan calados negativos (cuando el dominio superficial es menor al dominio de toda la malla).   
 !'eval' se ha evaluado inicialmente en la subrutinas aguassomerassubt
 !>Se aplica la condici�n: 
 !>Modelo superficial - velocidad nula en contorno m�vil 
 !'eval' s�lo se ha sobreescrito en la subrutina mallaaguassomeras para it=1.
 !Si alg�n coeficiente no es nulo se entra aqu�.
 !Si todos sus coeficientes son nulos (se ha considerado todo el dominio) y en alguna iteraci�n hay valores negativos de calado negativo,
 !se entrar� aqu� y 'eval' se modificar� dentro de la subrutina. 
 !Una vez que se entra siempre se entrar� a no ser que en nuevamalla se den valores nulos a todos los coeficientes.
 !>Modelo conjunto - no se aplica velocidad (posteriormente, en la subrutina mallaaguassomeras, se dan velocidades subterr�neas) en contorno m�vil 
 !Para it=1 se habr�n aplicado condiciones de velocidad nula (depende de la condici�n inicial o aproximaci�n inicial) en el contorno m�vil
 !En it=2 � it>2 del modelo conjunto siempre se habr�n aplicado valores de velocidad subterr�nea calculados con el modelo subterr�neo.
 !Aunque la malla sea mayor o menor se habr�n aplicado en todos los nodos del contorno m�vil al hacer 1 it por cada convergencia del modelo subterr�neo. 
 !'eval' se sobreescribe en la subrutina mallaaguassomeras para cada iteraci�n.
 !Si alg�n coeficiente no es nulo se entra aqu�.
 !Si todos sus coeficientes son nulos (se ha considerado todo el dominio) y hay valores negativos de calado negativo,
 !eval se modificar� en la subrutina mallaaguassomeras.		  
 call nuevamalla(modelo,i,j,z,vt,vth,vb,eval,vo)
 !As�, 'vo' seguramente ser� diferente la proxima vez que se resuelva una iteraci�n.		
endif
endif
deallocate(vn,vv,vec,voc,vu,vef,fx,fy,fz,vth,veac,veoc,veuc,ru,rb)	  
end

!----------------------------------------------------------------------------------------------------------------------------------------------------------
!Subrutina AGUASSUBTERRANEA (1 ecuaci�n de continuidad con las caracter�sticas de conductividad, �ngulo de anisotrop�a y porosidad). 
!En esta subrutina se calcula el nivel fre�tico (altura de la l�mina de agua). 
!Se contruye el sistema con la matriz de rigidez (una caja). Se calculan los niveles fre�ticos en los nodos esquina. 
!Las velocidades de Darcy se calculan postproceso. 
!Como soluci�n en cada iteraci�n habr� altura subterr�nea y una velocidad (en vtv=vbv y velx,vely).
!Se har�n todas las iteraciones no lineales para la convergencia del modelo.
!'it' es la iteraci�n del m�todo para resolver la no linealidad del sistema (s�lo se usa Picard). 
!'tol' gestiona la parada del bucle cuando dos soluciones consecutivas son similares.
!'vnv,vov,vuv' son contadores de nodos cuadr�ticos, lineales o con condici�n de contorno
!'veic' es el vector que contiene las integrales de contorno, 'vefv' es el t�rmino fuente de caudal por lluvia para flujo subterr�neo
!'qxx,qyy' llevan los valores de condiciones de contorno de caudal en los nodos donde �stas se conocen (guardadas en qx,qy), 
!y valores de caudales por metro lineal calculados a trav�s de las velocidades subterr�neas en el resto de nodos.
!'vecv,vocv' forman el t�rmino independiente del sistema  
!'vvv,vtv,vbv' son los vectores donde se guarda la soluci�n del sistema en cada iteraci�n, 'c' es el tama�o del sistema en cada iteraci�n.
!'inc' es la dimensi�n de la matriz del sistema cuadrada si se usasen elementos cuadr�ticos para niveles fre�ticos (n�ecuaciones*i)
!Otras variables con el mismo nombre son comentadas en la subrutina aguassomeras.                                                                              
!---------------------------------------------------------------------------------------------------------------------------------------------------------
subroutine aguassubterranea(i,j,x,y,z,zp,zzp,nd,kix,kiy,ag,sino,bcg,At,vnin,vov,eval,vb,ql,qb,qx,qy,vibv,vtv,velx,vely) 	 !*** A�adidos z,zzp,vb
use allocatacion
integer*4, dimension(:),allocatable::vnv
integer*4 i,j,u,uu,c,it,vnin(i),inc,ndim,k,conv			 
real*8, dimension(:),allocatable::vbdin,vecdin,vvv,vefv,vuv,vecv,vocv,vbv,qxx,qyy,veic		  	 
real*8 x(i),y(i),z(i),zp(i),zzp(i),tol,eval(i),At,vov(i),ag(i),alfa,vb(3*i)             !*** A�adidos z(i),zzp(i),alfa,vb(3*i)
real*8 vtv(i),vibv(i),kix(i),kiy(i),nd(i),qx(i),qy(i),qb(i),ql(i),velx(i),vely(i)									  		
logical nonzero	  
character sino*2,bcg*2

allocate(vnv(i),vvv(i),vefv(i),vuv(i),vecv(i),vocv(i),vbv(i),qxx(i),qyy(i),veic(i)) 

30  format(2(3x,A2,'(',I5,')=',E15.8E2))	  
	 
!Inizializaci�n previa de variables:
!-----------------------------------
inc=i
do u=1,i
vbv(u)=0.0_8
enddo
!*** Al subir alfa pasa mucho m�s caudal en la interfaz y tarda m�s iteraciones (probado de 0.00005 a 0.0005 - entre 0.0001 y 0.0005 hay mucho cambio y 
!*** para 0.001 dar�a error). En general alfa debe ser mayor que la conductividad m�nima y tanto mayor sea se debe trabajar con At menores al variar mucho 
!*** m�s la solucion.
!*** Para un gran paso de caudal en el embalse la l�mina toma curvatura para los n�meros de manning existentes
alfa=0.000025 !*** A�adida esta l�nea

!Inicio de iteraciones para un incremento de tiempo o para el problema estacionario (se aplica un esquema impl�cito):
!-------------------------------------------------------------------------------------------------------------------- 
it=0
tol=1e-6
do while (tol.eq.1e-6) 

!Selecci�n de en qu� iteraciones se utilizar� la soluci�n anterior para el m�todo de los gradientes biconjungados
!----------------------------------------------------------------------------------------------------------------
if (it.eq.0)then
nonzero=.false.	
else
nonzero=.true.
endif

!Inizializaci�n de variables (procedimientos ya aplicados con la subrutina aguassomeras):
!----------------------------------------------------------------------------------------
 do u=1,i
 vvv(u)=vtv(u)
 vocv(u)=0.0_8
 vefv(u)=0.0_8
 veic(u)=0.0_8
 enddo
 !El vn original es ahora vnin. El nuevo vn ser� diferente si el dominio es menor al de toda la malla.
 !En la ecuaci�n se eliminar�n las ecuaciones relativas a todos los nodos cuadr�ticos y no s�lo los del dominio suterr�neo.
 !Con vn se eliminar�n los de este dominio y el resto se eliminar�n con las CC al no considerar el resto de la malla.
 do u=1,i
  if ((vnin(u).eq.1).and.(eval(u).eq.0.0)) then
  vnv(u)=1
  else
  vnv(u)=0
  endif
 enddo
 !Se dimensionan ita y sa (tambi�n cia y ca), de acuerdo al n�mero total de coeficientes que se generar�n con las matrices elementales.
 !Adem�s ita ya lleva referenciados el n�mero de coeficientes que habr� por fila (en sus primeras i+1 componentes) antes de dar el formato MSR.
 !Cuando se calcule una matriz se almacenar�n sus coeficientes directamente en estos vectores.
 call dimvectsb (i,j,vnin)
 ndim=ita(inc+1)-1
	
!En caso de resolver de forma transitoria con incrementos de tiempo (esquema impl�cito considerando las matrices de masa).
!-------------------------------------------------------------------------------------------------------------------------
!Si se aplica el modelo conjunto y tiene afecci�n la condici�n similar a la condici�n seco-mojado, se debe tener en cuenta que las condiciones
!de contorno que se impondr�n en el contorno m�vil (que no coincidir�n con los valores de la soluci�n temporal previa y que son 
!generadas por el modelo superficial) se corresponden con un tiempo igual al tiempo para el que se busca la soluci�n.
if (sino.eq.'si') then
call timesubt (i,j,x,y,nd,At)
!Otra opci�n es utilizar la siguiente subrutina en vez de la anterior, calculando la matriz de masa concentrada:
!call timesubtconc(i,j,x,y,nd,At,ndim) !c�
do u=1,ndim
cia(u)=ita(u)
ca(u)=sa(u)
enddo
call orden (inc,k)
do u=1,i
 do uu=ita(u),ita(u+1)-1
 vocv(u)=vocv(u)+sa(uu)*vibv(ita(uu))
 enddo
 vocv(u)=vocv(u)+sa(u)*vibv(u)
enddo
do u=1,ndim
ita(u)=cia(u)
sa(u)=ca(u)
enddo
endif
do u=1,i
vtv(u)=0.0_8
enddo						

!Arreglo de las condiciones de contorno.
!---------------------------------------
!No es necesario sustituir las CC en el vector inicial ya que las integrales de contorno no dependen de las variables del sistema.
!Para el c�lculo de las integrales de contorno se usan los valores de caudal qx,qy (en los nodos lineales) le�dos inicialmente del fichero. 
!En vez de dar valores nulos donde no hay condici�n (por ejemplo poniendo qxx=0, qyy=0 si qx(u),qy(u)=sqrt(3)), se calculan unos valores aproximados de 
!caudal (con valores de velocidades calculados post-proceso en la iteraci�n anterior). No se tendr�an a trav�s de la soluci�n ya que esta variable no es 
!soluci�n del sistema. Sin embargo, �stos no tendr�n influencia alguna, pues los coeficientes calculados con ellos desaparecen por ensamblar 
!vectores elementales o por aplicar condiciones de nivel fre�tico sobre el sistema (en el contorno m�vil siempre se aplican). 	 
do u=1,i		   
 if (qx(u).eq.sqrt(3.0_8))then	   
  !Se tiene en cuenta la condici�n de espesor m�nimo que tambi�n aplica en la subrutina cajasasubt.
  !De este modo se evita utilizar espesores negativos cuando aparecen (teniendo esta condici�n un fin parecido al que se persigue con la aplicaci�n
  !de la condici�n seco-mojado para el modelo superficial).				
  if ((vvv(u)-zp(u)).le.0.0) then
  qxx(u)=velx(u)*0.001_8
  else
  qxx(u)=velx(u)*(vvv(u)-zp(u))
  endif
 else
 qxx(u)=qx(u)
 endif
 if (qy(u).eq.sqrt(3.0_8))then
  if ((vvv(u)-zp(u)).le.0.0) then
  qyy(u)=vely(u)*0.001_8
  else
  qyy(u)=vely(u)*(vvv(u)-zp(u))
  endif
 else
 qyy(u)=qy(u)
 endif
enddo

!C�lculo de otras matrices y vectores para la ecuaciones subterr�nea:
!--------------------------------------------------------------------
!En la siguiente subrutina se calculan las integrales de contorno.
call vectorcontornocaudales(i,j,x,y,zp,zzp,veic,qxx,qyy,vb,alfa)    !*** A�adidos zp,zzp,vb,alfa
!C�lculo de las integrales correspondientes a la lluvia.
call vectorcontornofuentesub(i,j,x,y,vefv,ql)
do u=1,i
vecv(u)=vocv(u)-veic(u)-qb(u)+vefv(u)
enddo
!C�lculo de la matriz A (Asx+Asy).
call cajasasubt (i,j,x,y,zp,vvv,kix,kiy,ag)

!Se inicializa vuv, para imponer las condiciones de contorno sobre el sistema
!----------------------------------------------------------------------------
do u=1,i
vuv(u)=vov(u)
enddo

!Aplicaci�n del m�todo de Picard (procedimientos ya aplicados con la subrutina aguassomeras):
!--------------------------------------------------------------------------------------------
call reducciondelsistema(i,c,vecv,vuv,vnv,vbv,inc,bcg,ndim)				
!Resolucion mediante gradientes biconjugados
allocate(vbdin(c),vecdin(c))  
do u=1,c			
vecdin(u)=vecv(u)
!Si (nonzero), que se da para it>0, vbv tendr� coeficientes no nulos. En otro caso siempre ser� nulo. 
vbdin(u)=vbv(u)
enddo
!Aqu� s�lo se usa el esquema impl�cito.
conv=1
u=0 
!Improbable que haya alg�n problema (conv ser� 0 tras llamar a cualquiera de las opciones para resolver el sistema)
if (bcg.eq.'no') then
 do while (conv.eq.1)
 call gradientesbiconjugados(vecdin,c,vbdin,nonzero,conv)
 nonzero=.true.
 u=u+1
   if ((u.eq.3).and.(conv.eq.1)) then
   write(6,*) 'No es posible obtener solucion del sistema en iteracion no lineal',it
   write(6,*) 'Pulsa ENTER. Si hay problemas a partir de aqui escoge otro metodo'
   read(5,*)
   exit
   endif
 enddo
else
 do while (conv.eq.1)
 call dslubc (vecdin,c,vbdin,ndim+7*c,ndim+3*c+12,conv)
 u=u+1
   if ((u.eq.3).and.(conv.eq.1)) then
   write(6,*) 'No es posible obtener solucion del sistema en iteracion no lineal',it
   write(6,*) 'Pulsa ENTER. Si hay problemas a partir de aqui escoge otro metodo'
   read(5,*)
   exit
   endif
 enddo		  
endif
do u=1,c			
vtv(u)=vbdin(u) 
enddo

!Arreglo del vector soluci�n introduciendo los valores prescritos de hd (y hd nulo en zonas donde hay soluci�n aguas someras
!tanto en nodos cuadr�ticos como lineales) y los valores nulos en los nodos cuadr�ticos (de la nueva malla).
!En caso de aplicar el modelo conjunto se tendr�n valores nulos de las variables donde no se calcule flujo subterr�neo.... si tiene afecci�n la cond
do u=1,i
 if (vuv(u).ne.sqrt(2.0_8)) then
  if (u.le.c) then
  do uu=c,u,-1
  vtv(uu+1)=vtv(uu)		   
  enddo
  endif
  c=c+1
  vtv(u)=vuv(u)
 endif 
enddo

!Postproceso tras la resoluci�n del sistema:
!-------------------------------------------
!Escritura de la soluci�n sobre el vector vbv.
do u=1,i
vbv(u)=vtv(u)
enddo
!C�lculo de las velocidades tras cada iteraci�n.
call velocidadessubterraneas (i,j,x,y,kix,kiy,ag,vtv,velx,vely)
it=it+1

!Se muestra la soluci�n de velocidades subterr�neo y niveles fre�ticos por pantalla (s�lo en los nodos del dominio subterr�neo si la 
!condici�n similar a la condici�n seco-mojado tiene afecci�n).	Impresi�n de resultados fraccionada, de 200 en 200
!l�neas (necesario definir 'ua' como entero).
!write(6,*)' '
!ua=1
!uu=0
!do u=1,i
! if ((vnin(u).ne.1).and.(eval(u).eq.0.0_8)) then		                     
! write(6,30)'us',u,velx(u),'vs',u,vely(u)
! uu=uu+1
!  if (uu.eq.200*ua) then
!  write(6,*)'Mostrados arriba 200 resultados. Pulsa enter.'
!  read(5,*)
!  ua=ua+1
!  endif
! endif	
!enddo
!write(6,*)' '
!ua=1
!uu=0
!do u=1,i
! if ((vnin(u).ne.1).and.(eval(u).eq.0.0_8)) then							 
! write(6,30)'hd',u,vtv(u)
! uu=uu+1
!  if (uu.eq.200*ua) then
!  write(6,*)'Mostrados arriba 200 resultados. Pulsa enter.'
!  read(5,*)
!  ua=ua+1
!  endif
! endif
!enddo

!Control para la parada del bucle para la iteraci�n no lineal.
write(6,*)'Error entre iteraciones: ',maxval(abs(vvv-vtv))
!Cuando la diferencia entre cualquier componente del vector soluci�n y la obtenida en la iteraci�n anterior sea menor que la variable tol
!se detendr�n las iteraciones (se sobreescribe aqu� tol y as� se saldr� del bucle dentro de esta subrutina).
if (maxval(abs(vvv-vtv)).lt.tol) then
tol=1.0
write(6,*)'Convergencia modelo subterraneo en iteracion:', it
endif

if (bcg.eq.'no') then
deallocate(vbdin,vecdin,ita,sa)
else
deallocate(vbdin,vecdin,cia,ca,cja)  
endif 
enddo

!Si se aplica el modelo conjunto y tuvo afecci�n la condici�n similar a la seco-mojado se habr�n dado condiciones de contorno
!en el contorno m�vil (previamente, a trav�s de la subrutina mallasubterranea, se dan niveles fre�ticos).
!As�, 'vov' seguramente ser� diferente la proxima vez que se resuelva esta ecuaci�n. 

!*** Conservaci�n de masa en el contorno m�vil
call conservacion (i,j,x,y,z,zp,zzp,vtv,velx,vely,vb) !*** A�adida esta l�nea
								  
deallocate(vnv,vvv,vefv,vuv,vecv,vocv,vbv,qxx,qyy,veic)
end

!----------------------------------------------------------------------------------------------------------------------------------------------
!Subrutinas de gesti�n del movimiento de contornos m�viles.
!----------------------------------------------------------------------------------------------------------------------------------------------
!Subrutina MALLAAGUASSOMERAS
!Esta subrutina puede modificar la malla del dominio superficial (sobreescribiendo el fichero malla.txt) haciendo una selecci�n del dominio 
!superficial en base al dominio subterr�neo, localizando el mismo contorno m�vil definido con la subrutina mallasubterranea.
!Adem�s da CC de velocidad subterr�nea (o nula si no la hay en la CI o aproximaci�n inicial o no se ha calculado soluci�n subterr�nea) 
!en el contorno m�vil.
!Se aplica antes de resolver/aplicar la ec superficial.
!Se utiliza una vez (1 it) si se aplica el modelo superficial para definir inicialmente el dominio superficial 
!dando CC de valor nulo. Las siguientes selecciones (por iteraci�n) del dominio superficial corresponden a la subrutina nuevamalla.
!Se utiliza en cada iteraci�n (en cada iteraci�n conjunta) si se aplica el modelo cunjunto (en la primera iteraci�n tambi�n se define el 
!dominio y se dan CC de valor nulo). 
!----------------------------------------------------------------------------------------------------------------------------------------------
subroutine mallaaguassomeras(i,j,z,zzp,vta,velx,vely,eval,vuc,vt)  !*** A�adido zzp
use interaccion
integer*4 i,j,u,ui,uuu,uu,io		   
real*8 z(i),zzp(i),vt(3*i),vta(i),velx(i),vely(i),eval(i),vuc(3*i),fi	 !*** A�adido zzp(i)											   
character fe*1

10  format(I5)
21  format(3/,I5)								
23  format(4/,A80)									
26  format(6X,6(X,I5))	
27  format(I5,X,6(X,I5))	
46  format(X,I5,X,A1,X,F11.7)						

!Inizializaci�n previa de variables:
!-----------------------------------
io=0
uu=0							  
do u=1,i
eval(u)=sqrt(3.0_8)
enddo
do u=1,3*i		  
vuc(u)=sqrt(2.0_8)	    				   
enddo

!Selecci�n de dominio considerando los valores de nivel fre�tico:
!----------------------------------------------------------------
open(unit=4,file='C:\mallainicial.txt',status='old')
!La opci�n scratch crea un archivo temporal que desaparece al cerrar el archivo
open(2,status= 'scratch') 
read(4,21)j
read(4,'(A)')a
 do u=1,j
 read(4,26) no(1),no(4),no(2),no(5),no(3),no(6)
   !Selecci�n inicial. Aqu� s�lo se consideran los elementos con todos sus nodos de nivel fre�tico nulo.   
   if ((vta(no(1)).eq.0.0_8).and.(vta(no(2)).eq.0.0_8).and.(vta(no(3)).eq.0.0_8)) then
   uu=uu+1
   write(2,27)uu,no(1),no(4),no(2),no(5),no(3),no(6)
   !Se considerar� el elemento.
   do uuu=1,6
   eval(no(uuu))=0.0_8
   enddo
   endif
 enddo
 rewind(4)
 read(4,23)a
 do u=1,j
 read(4,26)no(1),no(4),no(2),no(5),no(3),no(6)
   do ui=1,3
   b=ui+1-sb(ui)
   c=ui+2-sb(ui+1)
   !Se toman los elementos que rodean a la selecci�n inicial. Tambi�n se toman elementos pegados a nodos aislados o a l�neas que no forman elementos 
   !(regueros) donde no se calcul� soluci�n subterr�nea (niveles fre�ticos nulos). Es probable encontrar regueros pegados a la selecci�n inicial. 
   !Estos casos no se dar�n debido a una previa b�squeda de flujo superficial aislado en la soluci�n subterr�nea, ya que esto no se hace. 
   !Dentro de los elementos que se toman, no se ha calculado el flujo subterr�neo (hay soluci�n para el flujo superficial, de la anterior 
   !iteraci�n de la ecuaci�n superficial).   	   
   if ((vta(no(ui)).ne.0.0_8).and.(vta(no(b)).eq.0.0_8).and.(vta(no(c)).eq.0.0_8)) then  
	!Elementos con un nodo apoyado en el contorno m�vil - tienen nodos esquina con: dos niveles fre�tico nulos y uno no nulo.
	uu=uu+1
    write(2,27)uu,no(1),no(4),no(2),no(5),no(3),no(6)
	do uuu=1,6
    eval(no(uuu))=0.0_8
    enddo		 	
   elseif ((vta(no(ui)).ne.0.0_8).and.(vta(no(b)).ne.0.0_8).and.(vta(no(c)).eq.0.0)) then
	!Elementos con dos nodos apoyados en el contorno m�vil - tienen nodos esquina con: un nivel fre�tico nulo y dos no nulos. 
	!A trav�s de ellos se dan las CC de velocidad subterr�nea. No hay que considerar elementos con dos niveles fre�ticos nulos y uno no nulo 
	!para ello pues a�n en caso de que haya dos o m�s elementos juntos de este tipo compartir�n este nodo de nivel fre�tico no nulo, y este nodo 
	!formar� finalmente parte de un elemento de los analizados aqu�.
	uu=uu+1
    write(2,27)uu,no(1),no(4),no(2),no(5),no(3),no(6)
	do uuu=1,6
    eval(no(uuu))=0.0_8
    enddo
	 !Se llega a una interfaz sumando elementos (considerando tambi�n el caso particular, se hace a continuaci�n) a la selecci�n inicial y ser� la 
	 !interfaz que se calcul� en la subrutina mallasubterr�nea. Por tanto, la interfaz es la misma para flujo de aguassomeras y flujo subterr�neo. 
	 !Est� situada en nodos donde hay calado, en la orilla del dominio superficial. Por tanto se pasan valores que ya tienen esos nodos.
	 !Ser�n CC de velocidad subterr�nea. Se da un valor interpolado en los nodos cuadr�ticos (en ellos no se ha calculado valor).
	 if (vuc(no(ui+3)).ne.sqrt(2.0_8))then  
	 !Si se tienen dos elementos no particulares cuyo lado compartido tiene nodos lineales que se apoyan en contorno m�vil pero el cuadr�tico
	 !no lo hace, se dar�an condiciones de contorno con ambos elementos en dicho nodo. No se deben dar y aqu� se elimina la condici�n de contorno 
	 !la segunda vez que se vaya a dar. Si hay condici�n de contorno se elimina. Ning�n elemento particular que la elimine ser� considerado despu�s.
	 vuc(no(ui+3))=sqrt(2.0_8)		  
     vuc(i+no(ui+3))=sqrt(2.0_8)	      
	 else
	 vuc(no(ui))=velx(no(ui))   		  	 
     vuc(no(b))=velx(no(b))
	 vuc(no(ui+3))=0.5_8*(velx(no(ui))+velx(no(b)))
	 vuc(i+no(ui))=vely(no(ui))   		  	 
     vuc(i+no(b))=vely(no(b))
	 vuc(i+no(ui+3))=0.5_8*(vely(no(ui))+vely(no(b)))
	 endif
	 !Dado que la interfaz (contorno m�vil) es la misma, es posible garantizar la conservaci�n. Las velocidades subterr�neas estar�n mejor 
	 !aproximadas tanto menor sean los elementos pegados a la interfaz, dado que se interpolan linealmente.
	 !En otro caso (dos contornos m�viles) el n�mero de nodos en cada interfaz puede ser diferente, y aunque sea el mismo ser� complicado 
	 !calcular las componentes en nodos diferentes de forma que el caudal se conserve.	   
   endif
   enddo    
 enddo
 
!Ojo, es posible que el flujo subterr�neo haga casi desaparecer al flujo superficial, si la dimensi�n de �ste es peque�a (por ejemplo 
!un r�o de poco de ancho). En este caso, podr�a seleccionarse el dominio a partir de un reguero, y este reguero estar dividido en dos 
!partes en alguna zona. En este caso habr� que refinar m�s la malla.
!A�n en caso de que los elementos considerados pegados a ambos trozos de reguero se toquen formando un dominio cont�nuo habr�a problemas
!ya que se dar� condici�n de velocidad nula en el punto donde ambas partes se unen.

!Caso particular:
!----------------
!Es posible tener elementos que no se consideran con el procedimiento anterior y que son necesarios para tomar todos los elementos hasta
!el contorno m�vil.
!Estos elementos no estar�n pegados a la selecci�n inicial, y cumplen la particularidad de tener sus tres nodos esquina apoyados en el contorno m�vil.
!Podr�an tampoco estar pegados a la �ltima selecci�n hecha si forman una banda con el ancho de un elemento.  
rewind(4)
read(4,23)a
 do u=1,j
 read(4,26) no(1),no(4),no(2),no(5),no(3),no(6)   
   !Se selecciona el elemento (pegado o no pegado a la �ltima selecci�n hecha antes de entrar aqu�).
   if ((z(no(1)).eq.zzp(no(1))).and.(z(no(2)).eq.zzp(no(2))).and.(z(no(3)).eq.zzp(no(3)))) then	  !*** Cambiados los dos if que hab�a por este if
   uu=uu+1
   write(2,27)uu,no(1),no(4),no(2),no(5),no(3),no(6)		
   do uuu=1,6
   eval(no(uuu))=0.0_8
   enddo 
    !Se corrijen las condiciones de contorno. Se consideran diferentes casos en funci�n de si el elemento
	!est� o no pegado a la �ltima selecci�n hecha (selecci�n que puede estar continuamente cambiando).
	do ui=1,3
    b=ui+1-sb(ui)
    c=ui+2-sb(ui+1)
    e=ui+4-sb(ui)
    f=ui+5-sb(ui+1) 
	 !Elemento no pegado a la �ltima selecci�n (pegados por 2 nodos o 1 nodo sin lado pegado o totalmente separados).
	 if ((vuc(no(ui+3)).eq.sqrt(2.0_8)).and.(vuc(no(e)).eq.sqrt(2.0_8)).and.(vuc(no(f)).eq.sqrt(2.0_8))) then
	 vuc(no(ui))=velx(no(ui))
	 vuc(i+no(ui))=vely(no(ui))
	 vuc(no(b))=velx(no(b))
	 vuc(i+no(b))=vely(no(b))
	 vuc(no(c))=velx(no(c))
	 vuc(i+no(c))=vely(no(c))
	 vuc(no(ui+3))=0.5_8*(velx(no(ui))+velx(no(b)))
     vuc(i+no(ui+3))=0.5_8*(vely(no(ui))+vely(no(b)))
	 vuc(no(e))=0.5_8*(velx(no(b))+velx(no(c)))
     vuc(i+no(e))=0.5_8*(vely(no(b))+vely(no(c)))
	 vuc(no(f))=0.5_8*(velx(no(ui))+velx(no(c)))
     vuc(i+no(f))=0.5_8*(vely(no(ui))+vely(no(c)))
	 exit
	 !Elemento pegado a la �ltima selecci�n hecha (pegado por uno de sus lados).
	 elseif ((vuc(no(ui+3)).ne.sqrt(2.0_8)).and.(vuc(no(e)).eq.sqrt(2.0_8)).and.(vuc(no(f)).eq.sqrt(2.0_8))) then
	 vuc(no(c))=velx(no(c))
	 vuc(i+no(c))=vely(no(c))
	 vuc(no(ui+3))=sqrt(2.0_8)
     vuc(i+no(ui+3))=sqrt(2.0_8)
	 vuc(no(e))=0.5_8*(velx(no(b))+velx(no(c)))
     vuc(i+no(e))=0.5_8*(vely(no(b))+vely(no(c)))
	 vuc(no(f))=0.5_8*(velx(no(ui))+velx(no(c)))
     vuc(i+no(f))=0.5_8*(vely(no(ui))+vely(no(c)))
     exit
	 !Elemento encerrado por la �ltima selecci�n hecha (pegado por dos de sus lados).
	 elseif ((vuc(no(ui+3)).ne.sqrt(2.0_8)).and.(vuc(no(e)).ne.sqrt(2.0_8)).and.(vuc(no(f)).eq.sqrt(2.0_8))) then
	 vuc(no(ui+3))=sqrt(2.0_8)
     vuc(i+no(ui+3))=sqrt(2.0_8)
	 vuc(no(e))=sqrt(2.0_8)
     vuc(i+no(e))=sqrt(2.0_8)
	 vuc(no(f))=0.5_8*(velx(no(ui))+velx(no(c)))
     vuc(i+no(f))=0.5_8*(vely(no(ui))+vely(no(c)))
	 exit
	 endif
    enddo
   endif
   !*** Eliminado uno de los endif que hab�a
 enddo

!Valores de la iteraci�n anterior:
!---------------------------------
!Se toman como valores de la iteraci�n anterior (variables de la ecuaci�n superficial) los valores de las CC de velocidad dados.
!Se trata valores en nodos del contorno m�vil. De esta forma afectar�n a las integrales de contorno (y existir� conservaci�n).
!Ojo, esto da problemas en la convergencia (al modificar la variable 'vt') si se usa la primera opci�n propuesta a continuaci�n.
do u=1,2*i
 if (vuc(u).ne.sqrt(2.0_8)) then		   
 vt(u)=vuc(u)
 endif
enddo

!Arreglo de las condiciones de contorno:
!---------------------------------------
!Fuera de la nueva malla se dan condiciones de contorno nulas.
!El sistema tendr� coeficientes s�lo referidos los nodos de la malla que se selecciona, pero el sistema tendr� la dimensi�n correspondiente 
!a los nodos de toda la malla (no se pierde la referencia del nodo). Con esto tambi�n se reducir� el tama�o del sistema para considerar el 
!sistema correspondiente.
do u=1,i
if (eval(u).eq.sqrt(3.0_8)) then
vuc(u)=0.0_8
vuc(i+u)=0.0_8
vuc(2*i+u)=0.0_8
endif
enddo

!Introducci�n de las CC de flujo superficial que existen en el contorno no m�vil (le�das desde fichero). El contorno no m�vil ser� el trozo de 
!contorno que tambi�n encierra a la selecci�n de elementos y que es parte del contorno de toda la malla (del dominio conjunto).
do u=1,i+2
 read(4,'(A)')a
enddo
eaf=0
do while (eaf.eq.0)
 read (4,46,iostat=eaf) u,fe,fi
 io=io+1 
!c�
!Intersecci�n del contorno m�vil con el contorno del dominio conjunto (existe CC en el fichero). 
!>Primera opci�n que permite la utilizaci�n exacta de las CC del fichero. En la intersecci�n:
!Se deja el valor de velocidad subterr�nea �nicamente como valor de la iteraci�n anterior (se usan en las integrales de contorno). 
!Se da la condici�n de contorno existente (de velocidad o de calado, se usan en la reducci�n del sistema).
!Respecto a la convergencia, en la itersecci�n, si hay CC de calado la velocidad calculada se modifica evitando 
!la convergencia y si hay CC de velocidad la velocidad de la iteraci�n anterior nunca tiene ese valor. 
!El problema se evita modificando la variable vt debajo, fuera de este bucle "do while", y no arriba (si se modifican valores de la iteraci�n anterior, 
!vuc y vt deben coincidir para que haya convergencia). Pero no se utilizar�n las velocidades subterr�neas en la intersecci�n (ojo con la conservaci�n).
! if (eaf.eq.0) then 
!  if (io.le.i) then
!   !Hay CC de velocidad en fichero, se sobreescribe la condici�n en caso de que se haya dado condici�n de velocidad (intersecci�n).
!   !En otro caso vuc tambi�n se sobreescribir� (pues vale sqrt(2.0_8) indicando que no hay condici�n). 
!   if ((vuc(u).ne.0.0_8).and.(fe.ne.'o')) then 
! 	vuc(u)=fi
! 	endif
!  elseif ((io.gt.i).and.(io.le.2*i)) then
!   !Lo mismo pero para la velocidad en direcci�n y.
!   if ((vuc(i+u).ne.0.0).and.(fe.ne.'o')) then 
! 	vuc(i+u)=fi
! 	endif															   
!  else
! 	!Hay CC de calado en fichero, se elimina la condici�n en caso de que se haya dado condici�n de velocidad (intersecci�n).
! 	!Se da el valor sqrt(2.0_8) para indicar que no hay condici�n en todos los casos, incluyendo al caso anterior (donde se da CC de calado
!   !en el fichero, nunca se da CC de velocidad).
! 	if ((vuc(2*i+u).ne.0.0_8).and.(fe.ne.'o')) then 	 
! 	vuc(u)=sqrt(2.0_8)
! 	vuc(i+u)=sqrt(2.0_8)
! 	vuc(2*i+u)=fi+z(u)
! 	endif
!  endif
! endif 
!>Segunda opci�n que asegura la convergencia pero sin utilizar exactamente las CC del fichero. En la intersecci�n:
!Si existe CC calado se a�ade CC de velocidad subterr�nea (dar s�lo la CC velocidad es complicado porque habr�a que dar CC en nodos cuadr�ticos 
!no pertenecientes al contorno m�vil).
!Si existe CC de velocidad se substituye por la de velocidad subterr�nea.
 if (eaf.eq.0) then 
  !Si vuc es distinto de cero (o eval es distinto de sqrt(3.0)) el nodo pertenece a la nueva malla.
  !Si vuc es sqrt(2.0_8) el nodo no tiene valor de contorno (nunca lo tendr� si no se est� en el contorno). 
  if (io.le.i) then
    !Hay CC de velocidad en el fichero. Se asegura que s�lo exista la CC dada de velocidad.
    if ((vuc(u).eq.sqrt(2.0_8)).and.(fe.ne.'o')) then 
 	vuc(u)=fi
 	endif
  elseif ((io.gt.i).and.(io.le.2*i)) then
    !Lo mismo pero para la velocidad en direcci�n y.
    if ((vuc(i+u).eq.sqrt(2.0_8)).and.(fe.ne.'o')) then 
 	vuc(i+u)=fi
 	endif															   
  else
    !Ocurre en todo el contorno m�vil que vuc(2*u+i) es sqrt(2.0_8) al no darse CC de calado a trav�s de la soluci�n subterr�nea.
 	!Hay CC de calado en el fichero. Se da la CC superficial en la intersecci�n y no se elimina la que se haya dado de velocidad
	if ((vuc(2*i+u).eq.sqrt(2.0_8)).and.(fe.ne.'o')) then 	
 	vuc(2*i+u)=fi+z(u)
 	endif
  endif
 endif
enddo
		
!Escritura de la malla:
!----------------------
!La malla tendr� la i original porque los n�meros de nodo ir�n hasta i, y la nueva j para construir bien las cajas existentes
!(no se pierde la referencia del nodo).
rewind(2)                                 
close(4)       
open(1,file = 'C:\malla.txt',status='replace')
write(1,'(A)')'nudos'
write(1,10)i
write(1,'(A)')'elementos'
write(1,10)uu
write(1,'(A)')'elemento p1  p2  p3  p4  p5  p6 (antihorarios)'
do u=1,uu
 read(2,'(A)')a
 write(1,'(A)')a
enddo
close(2) 
close(1)
j=uu
end

!----------------------------------------------------------------------------------------------------------------------------------------------
!Subroutina NUEVAMALLA
!Esta subrutina aplica la condici�n seco-mojado o la condici�n similar a la condici�n seco-mojado. 
!Se aplica tras resolver la ec superficial.
!>Se utiliza (condici�n seco-mojado) en cada iteraci�n si se aplica el modelo superficial para definir el dominio superficial. 
!Se hace la selecci�n o deselecci�n de elementos para el dominio superficial (en malla.txt) en funci�n de los calados calculados con la 
!ecuaci�n correspondiente y se dan nuevos valores de calado como valores de la iteraci�n anterior. Se da CC de valor nulo en el contorno m�vil.  
!>Se utiliza en cada iteraci�n si se aplica el modelo conjunto pero s�lo aplica valores de calado de la iteraci�n anterior para una posterior 
!selecci�n del dominio subterr�neo con la subrutina mallasubterr�nea. Cualquier selecci�n de elementos no se considera (no va en malla.txt). 
!Posibilidades (en esta iteraci�n, esto es, cada vez que se entra):
!Selecci�n de elementos pegados al contorno m�vil actual (u original). 
!Deselecci�n de cualquier elemento dentro del dominio delimitado por el contorno m�vil actual.
!Si el calado es positivo �nicamente en los nodos de una l�nea (regero), �stos no se consideran (habr�a que refinar la malla).
!----------------------------------------------------------------------------------------------------------------------------------------------
subroutine nuevamalla(modelo,i,j,z,vt,vv,vb,eval,vuc)
use interaccion
integer*4, dimension(:),allocatable::pa,pe
integer*4 i,j,u,ui,uuu,uu,io
real*8, dimension(:),allocatable::si								  
real*8 vt(3*i),vv(i),vb(3*i),eval(i),z(i),vuc(3*i),fi,hmin 
character modelo*12,fe*1

allocate(si(3*i),pa(i),pe(i))

10  format(I5)
21  format(3/,I5)
23  format(4/,A80)
26  format(6X,6(X,I5))	
27  format(I5,X,6(X,I5))
46  format(X,I5,X,A1,X,F11.7)
 												   
!Par�metro para evitar inestabilidades:
!--------------------------------------
!Se definen los elementos interfaz como aquellos elementos pegados al contorno m�vil y no considerados en el dominio.
!Es posible que la nueva soluci�n conlleve a deseleccionar elementos seleccionados en una iteraci�n anterior. 
!Este proceso puede ser indefinido si la altura del agua coincide aproximadamente con la cota de los nodos de los elementos interfaz. 
!La soluci�n es usar un parametro 'hmin' que asegure que exista una cierta altura para la selecci�n (se limita la adicci�n de elemetos).   
!El parametro deber� ser mayor tanto menor sea el n�mero de elementos interfaz.	
hmin=1e-3   
!hmin=2.0_8	!es

!Inizializaci�n previa de variables:
!-----------------------------------
io=0
uu=0
do u=1,3*i		  
vuc(u)=sqrt(2.0_8)
si(u)=0.0_8 
enddo
do u=1,i
eval(u)=sqrt(3.0_8)
enddo

!Evaluci�n de la soluci�n obtenida:
!---------------------------------- 			
!Los nuevos valores dados como valores de la iteraci�n anterior siempre ser�n err�neos (son supuestos) y m�s en velocidades al ser m�s variables. 
!Sin embargo es necesario dar valores de la iteraci�n anterior de calado (y altura) en los nuevos nodos dado que al aplicar la ecuaci�n de aguas 
!someras no puede haber valores nulos o negativos. Se dan en nuevo contorno m�vil si se a�aden elementos (hay valores nulos), y no se dan en �l
!si se eliminan elementos (hay valor). 
!No es necesario dar valores de la iteraci�n anterior de velocidad excepto en el nuevo contorno m�vil donde tomar�n el mismo valor que las CC.
!Por tanto no se dar�n en el contorno m�vil actual ni en nodos cuadr�ticos entre contorno m�vil actual y nuevo si se a�aden elementos (hay valores 
!nulos o de velocidad subterr�nea en el primer caso y valores nulos en el segundo caso), ni se dan en cualquier nodo perteneciente al dominio pero 
!no al nuevo contorno m�vil si se eliminan elementos (hay valor). 
open(unit=4,file='C:\mallainicial.txt',status='old') 
open(2,status= 'scratch')
read(4,21)j
read(4,'(A)')a

!Se decide si con el calado calculado se utiliza el siguiente elemento, no se utiliza, o se eliminan elementos de la malla anterior.
!Esto permitir�, por ejemplo, a�adir elementos en una parte del contorno, eliminar elementos en otra parte, y no hacer nada en otra parte. 
!Se trabaja con vv que es vt(2*i+u) original ya que vt(2*i+u) se sobreescribir�. 
  do u=1,j
  read(4,26)no(1),no(4),no(2),no(5),no(3),no(6)
  !Primer proceso. Aqu� s�lo se consideran los elementos con todos sus nodos positivos. O se seleccionan los mismos elementos de la malla anterior 
  !o se seleccionan menos elementos (eliminaci�n de elementos en este caso). 
   if ((vv(no(1)).gt.0.0_8).and.(vv(no(2)).gt.0.0_8).and.(vv(no(3)).gt.0.0_8)) then  
    uu=uu+1
    write(2,27)uu,no(1),no(4),no(2),no(5),no(3),no(6)
	cycle 
   endif 
   !Segundo proceso. Se analizan elementos pegados al contorno m�vil original que no pertenecen a la malla anterior. O no se selecciona o se 
   !selecciona (adicci�n en este caso). 
   do ui=1,3
    b=ui+1-sb(ui)
    c=ui+2-sb(ui+1)
    e=ui+4-sb(ui)
    f=ui+5-sb(ui+1)
	if ((vv(no(ui)).gt.0.0).and.(vv(no(b)).gt.0.0).and.(vv(no(c)).eq.0.0)) then
	!Elementos con dos nodos apoyados en el contorno m�vil original. Si la altura en ambos es superior a la cota en el tercero (y 'hmin') se
    !considera el elemento.
    !Se dan valores como valores de la iteraci�n anterior. En el nuevo nodo se suman los valores de altura de los otros dos, a posibles valores ya 
	!tenidos en cuenta en este proceso. Por ello con la variable 'si' se realizar� la media de estas alturas en cada nodo posteriormente.
	if (((vv(no(ui))+z(no(ui))).gt.(z(no(c))+hmin)).and.((vv(no(b))+z(no(b))).gt.(z(no(c))+hmin))) then		   	  		 
	uu=uu+1
	write(2,27)uu,no(1),no(4),no(2),no(5),no(3),no(6)
	vb(2*i+no(c))=vb(2*i+no(c))+vv(no(ui))+z(no(ui))+vv(no(b))+z(no(b))	  
	si(2*i+no(c))=si(2*i+no(c))+2.0_8	
	exit
	endif
    !Elementos con un nodo apoyado en el contorno m�vil original. Si la altura en �l es superior a la cota en los otros dos nodos del elemento (y 'hmin')
	!se considera el elemento.
	!Se dan valores como valores de la iteraci�n anterior. En los dos nuevos nodos se suman el valor de altura del otro, a posibles valores ya tenidos 
	!en cuenta en este proceso. Por ello con la variable 'si' se realizar� la media en cada nodo posteriormente.   
    elseif ((vv(no(ui)).gt.0.0).and.(vv(no(b)).eq.0.0).and.(vv(no(c)).eq.0.0)) then
    if (((vv(no(ui))+z(no(ui))).gt.(z(no(b))+hmin)).and.((vv(no(ui))+z(no(ui))).gt.(z(no(c))+hmin))) then		 
	uu=uu+1
	write(2,27)uu,no(1),no(4),no(2),no(5),no(3),no(6)
	vb(2*i+no(b))=vb(2*i+no(b))+vv(no(ui))+z(no(ui))  
	vb(2*i+no(c))=vb(2*i+no(c))+vv(no(ui))+z(no(ui))  	 
	si(2*i+no(b))=si(2*i+no(b))+1.0_8
	si(2*i+no(c))=si(2*i+no(c))+1.0_8		
	exit
	endif						   	
	endif
   enddo
 enddo
 !C�lculo de valores medios:
 do u=1,i
 if (si(2*i+u).ne.0.0_8) then
 vb(2*i+u)=vb(2*i+u)/si(2*i+u)
 vt(2*i+u)=vb(2*i+u)-z(u)
 endif
 enddo
 !A continuaci�n se sobreescribe la variable 'eval' para que indique cu�les son los elementos considerados que forma la nueva malla, 
 !mayor, igual o menor a la anterior.
 rewind(2)
 do u=1,uu
 read(2,26)no(1),no(4),no(2),no(5),no(3),no(6)
   do uuu=1,6
   eval(no(uuu))=0.0_8
   enddo									    
 enddo

!Se aplica la condici�n seco-mojado (condici�n completa):
!--------------------------------------------------------
!Siempre necesaria cuando se hagan al menos 2 iteraciones consecutivas de la ecuaci�n de aguas someras (siempre en el modelo superficial). 
if (modelo.eq.'superficial') then  
 !Tercer proceso. Se seleccionan otros elementos que tambi�n pertenecer�an a la malla (hay calados positivos en los tres nodos esquina).
 rewind(4) 
 read(4,23)a
 do u=1,j
 read(4,26)no(1),no(4),no(2),no(5),no(3),no(6)
   if ((eval(no(1)).eq.0.0_8).and.(eval(no(2)).eq.0.0_8).and.(eval(no(3)).eq.0.0_8)) then	
	do ui=1,3
    b=ui+1-sb(ui)
    c=ui+2-sb(ui+1)
    e=ui+4-sb(ui)
    f=ui+5-sb(ui+1) 
    if (eval(no(ui+3)).ne.0.0_8) then
	uu=uu+1
	write(2,27)uu,no(1),no(4),no(2),no(5),no(3),no(6)
	!Se considera el elemento (s�lo necesario dar valores a 'eval' en nodos cuadr�ticos).
	eval(no(ui+3))=0.0_8
	eval(no(e))=0.0_8
	eval(no(f))=0.0_8	
	 !Se dar�n dos casos bajo esta condici�n:
	 !if ((si(2*i+no(ui)).ne.0).and.(si(2*i+no(b)).ne.0).and.(si(2*i+no(c)).ne.0)) then 	 
	  !Selecci�n de elementos encerrados.
	  !Al a�adir elementos pegados al contornos m�vil en zonas donde �ste tiene gran curvatura pueden quedar atrapados 
	  !elementos de este tipo entre el nuevo contorno m�vil. Se localizar� un contorno (con el algoritmo donde se da vuc=0) que los tiene en cuenta. 	 
	  !Si no se seleccionasen se trabajar�a con una malla de menos elementos donde algunos nodos del contorno no tendr�n condici�n.
	  !Proceso necesario aunque supondr� un error si hay alturas muy distintas en los nodos de cada elemento a a�adir.
	 !continue
	 !else
	  !Selecci�n de otros elementos particulares.
	  !Pueden aparecer varios elementos con un nodo esquina com�n fuera del contorno m�vil original de modo que s�lo se a�ade alguno de ellos. 
	  !Proceso no necesario, pues siempre se dar� condici�n en el contorno de la malla seleccionada.	   
	  !Se a�aden elementos con poco sentido f�sico cuando la discretizaci�n no es buena (poco refinamiento). Siempre supondr� 
	  !un error ya que la altura ser� muy distinta en los nodos de cada elemento a a�adir. Una posible soluci�n en este caso es subir 'hmin'.
	 !continue 
	 !endif
	exit
	endif
	enddo
   endif									    
 enddo

 !Se sobreescribe 'vuc' para dar condici�n de contorno de no deslizamiento (velocidades nulas).  
 rewind(4)
 read(4,23)a
 do u=1,j
 read(4,26)no(1),no(4),no(2),no(5),no(3),no(6)				 
   !Aqu� se localizan los nuevos elementos interfaz. Considerando que los calados pueden haber sido modificados, �stos ser�n  
   !Los elementos donde hay dos nodos de calado positivo y uno de calado negativo (uno o m�s elementos considerados en la iteraci�n anterior ahora 
   !se eliminan). �nico tipo de elementos si es la primera vez que se selecciona un dominio.
   !Los elementos donde haya dos nodos de calado positivos y uno de calado nulo (elemento ya considerado en la malla anterior o elemento a�adido).
   !El objetivo es localizar los nodos del contorno m�vil nuevo. Por ello la localizaci�n de elementos se hace mejor a trav�s de 'eval'. 
   !Con ello se asegura que se est�n tratando nodos del contorno m�vil nuevo. 
   do ui=1,3
   b=ui+1-sb(ui)
   c=ui+2-sb(ui+1)
   if ((eval(no(ui)).ne.sqrt(3.0_8)).and.(eval(no(b)).ne.sqrt(3.0_8)).and.(eval(no(c)).eq.sqrt(3.0_8))) then	 						 
	vuc(no(ui))=0.0_8   		  	 
    vuc(no(b))=0.0_8
	vuc(no(ui+3))=0.0_8
	vuc(i+no(ui))=0.0_8   		  	 
    vuc(i+no(b))=0.0_8
	vuc(i+no(ui+3))=0.0_8   	 	   
   endif
   enddo 
 enddo 

 !Cuarto proceso. Arreglo de condiciones de contorno en elementos casi despegados (pegados por un nodo) ya considerados. 
 !Es posible a�adir un elemento (o m�s) tal que el nuevo contorno m�vil haga una curva cerrada con �l (o ellos). 
 !Esto puede ocurrir para un (dos o m�s) elemento si hay al menos tres (cuatro o m�s) elementos fuera del dominio original que compartan un nodo del 
 !contorno m�vil y no tienen ning�n otro nodo en este contorno. Si se tienen tres elementos as� y se selecciona el elemento del medio (con los pasos
 !anteriores) aparecer� un problema. Con otra disposici�n geom�trica el tercer proceso a�adir� elementos a los lados de este elemento.  
 !El problema es que el elemento tendr� vuc=0 en todos sus nodos. As�, tanto el gradiente de la velocidad como la velocidad ser�n nulos y en la ecuaci�n 
 !de continuidad no se aportan coeficientes en las filas y columnas correspondientes a los nuevos dos nodos (al menos sin estabilizaci�n).   
 do u=1,i
 pa(u)=0
 pe(u)=0
 enddo
 rewind(4)
 read(4,23)a
 !Selecci�n de los nodos que causan el problema (para uno o m�s elementos de este tipo).
 do u=1,j
 read(4,26)no(1),no(4),no(2),no(5),no(3),no(6)				 
   if ((vuc(no(1)).eq.0.0_8).and.(vuc(no(2)).eq.0.0_8).and.(vuc(no(3)).eq.0.0_8)) then
   if ((vuc(no(4)).eq.0.0_8).and.(vuc(no(5)).eq.0.0_8).and.(vuc(no(6)).eq.0.0_8)) then
   !Si todos los nodos del elemento tienen condici�n de no delizamiento se sobrescribe la variable 'pa' en los nodos esquina.  
   !De momento la variable 'vuc' s�lo es nula en los nodos del contorno m�vil	 						 
	do ui=1,3
    pa(no(ui))=1
    enddo
	!Se evita que este elemento sea considerado como parte de la malla durante este proceso. As�, se trata una malla modificada.  	 	   
    cycle
   endif
   endif
   if ((eval(no(1)).eq.0.0_8).and.(eval(no(2)).eq.0.0_8).and.(eval(no(3)).eq.0.0_8)) then
   !Si el elemento est� en la malla modificada se sobreescribe la variable 'pe' en los nodos esquina.
   do ui=1,3
   pe(no(ui))=1
   enddo
   endif 
 enddo
 do u=1,i
 !El problema se resuelve dando condiciones de altura en los dos nuevos nodos citados. As� no existir�n esas ecuaciones elementales de continuidad.
 if ((pa(u).eq.1).and.(pe(u).eq.0)) then
 !Aqu� se seleccionan los nuevos nodos. Para entrar aqu� cada nodo debe pertenecer a un elemento con condici�n de contorno de no deslizamiento 
 !en todos los nodos y debe pertenece a un elemento que no est� en la malla modificada.
 !Un valor nulo de condici�n de altura en vez de el vb(2*i+u), supuesto como valor de la iteraci�n anterior, puede dar problemas de deselecci�n 
 !y selecci�n de estos elementos. 
 vuc(2*i+u)=vb(2*i+u) 		  	 
 endif
 enddo
 
 !En el nuevo contorno m�vil se dan valores nulos de velocidad. As�, se calcular�n bien las integrales de contorno.
 do u=1,2*i
  if (vuc(u).ne.sqrt(2.0_8)) then		   
  vt(u)=vuc(u)
  endif
 enddo

 !Fuera de la nueva malla considerada (en todo el dominio que no es de aguas someras) se dan condiciones de contorno nulas 
 !de velocidad y altura (u, v, Ht). As�, se reducir� el sistema eliminando filas y columnas que ya ser�n nulas (no se generar�n coeficientes para
 !elementos que no son de esta malla). Finalmente no se calcular� la soluci�n en esos nodos.   
 do u=1,i
  if (eval(u).eq.sqrt(3.0_8)) then
  vuc(u)=0.0_8
  vuc(i+u)=0.0_8
  vuc(2*i+u)=0.0_8
  endif
 enddo

 !Introducci�n de las CC de flujo superficial que existen en el contorno no m�vil (le�das desde fichero).
 do u=1,i+2
  read(4,'(A)')a
 enddo
 eaf=0
 do while (eaf.eq.0)
  read (4,46,iostat=eaf) u,fe,fi
  io=io+1  
  if (eaf.eq.0) then
   if (io.le.i) then
    if ((vuc(u).eq.sqrt(2.0_8)).and.(fe.ne.'o')) then
    vuc(u)=fi	   
    endif
   elseif ((io.gt.i).and.(io.le.2*i)) then
    if ((vuc(i+u).eq.sqrt(2.0_8)).and.(fe.ne.'o')) then
    vuc(i+u)=fi
    endif															   
   else
   !Al no darse condiciones de contorno de calado (vuc(2*i+u).eq.sqrt(2.0_8)) es equivalente a (eval(u).ne.sqrt(3.0)).
    if ((vuc(2*i+u).eq.sqrt(2.0_8)).and.(fe.ne.'o')) then	
    vuc(2*i+u)=fi+z(u)
    endif
   endif
  endif 
 enddo

 !Escritura de la nueva malla:  
 rewind(2)                             
 close(4)           
 open(1,file = 'C:\malla.txt',status='replace')
 write(1,'(A)')'nudos'
 write(1,10)i
 write(1,'(A)')'elementos'
 write(1,10)uu
 write(1,'(A)')'elemento p1  p2  p3  p4  p5  p6 (antihorarios)'
 do u=1,uu
  read(2,'(A)')a
  write(1,'(A)')a
 enddo
 close(2)
 close(1)
 j=uu

!Se aplica la condici�n similar a la condici�n seco mojado:
!----------------------------------------------------------
!Se aplica en otro caso (modelo conjunto). S�lo considera el primer y el segundo procesos.
!Con la subrutina mallaaguassomeras se llega por otro proceso a la selecci�n que produce el tercer proceso de esta subrutina.
!No es necesario arreglar 'vuc' ya que se van a calcular velocidades subterr�neas. �stas se introducen en la subrutina mallaaguassomeras.
!Al utilizar la subrutina mallaaguassomeras habr� velocidades subterr�neas en todos los nodos por lo que no ser� necesario el cuarto proceso
!de esta subrutina.
!Tampoco ser� necesario arreglar las condiciones de contorno o escribir la malla ya que se har� en la subrutina mallaaguassomeras.
else
 close(2)                             
 close(4)
 !La subrutina mallasubterranea se basar� para localizar el contorno m�vil (el mismo que localizar�a la condici�n seco-mojado) en los 
 !valores nulos o negativos de calado. Podr�a ocurrir que quedasen nodos con calado positivo fuera del nuevo contorno m�vil. 
 !Por ejemplo si en el primer proceso de esta subrutina no se considera un elemento que tiene dos nodos con calado positivo y uno de ellos en el 
 !nuevo contorno m�vil. Por tanto:
 do u=1,i
  if ((eval(u).eq.sqrt(3.0_8)).and.(vt(2*i+u).gt.0.0_8)) then
  vt(2*i+u)=0.0_8
  endif
 enddo
 !Obviamente, esto no es necesario si se calcula a continuaci�n otra iteraci�n de la ecuaci�n de aguas someras (con la condici�n seco-mojado).
endif
     
deallocate(si,pa,pe)
end

!---------------------------------------------------------------------------------------------------------------------------------------------
!Subrutina MALLASUBTERRANEA
!Esta subrutina puede modificar la malla del dominio subterr�neo (sobreescribiendo el fichero mallasub.txt) haciendo una selecci�n del dominio 
!subterr�neo en base al dominio superficial. 
!Adem�s da CC de altura, esto es, de nivel fre�tico en el contorno m�vil.
!Se aplica antes de resolver la ec subterr�nea y es equivalente a la subrutina mallaaguassomeras.
!S�lo se utiliza si se resuelve el modelo conjunto.
!Se utiliza una vez por cada grupo de iteraciones hasta convergencia del modelo subterr�neo (en cada iteraci�n conjunta) si se aplica el 
!modelo cunjunto. 
!---------------------------------------------------------------------------------------------------------------------------------------------
subroutine mallasubterranea (i,j,z,zp,zzp,vt,eval,qx,qy,vuc)
use interaccion
integer*4 i,j,u,ui,uuu,uu,io		   
real*8 vt(3*i),z(i),zp(i),zzp(i),eval(i),vuc(i),qx(i),qy(i),fi
character fe*1											 

10  format(I5)
21  format(3/,I5)								
23  format(4/,A80)									
26  format(6X,6(X,I5))	
27  format(I5,X,6(X,I5))	
46  format(X,I5,X,A1,X,F11.7)						
 
!Inizializaci�n previa de variables:
!-----------------------------------
io=0
uu=0						  
do u=1,i					 
qx(u)=sqrt(3.0_8)					 
qy(u)=sqrt(3.0_8)
eval(u)=sqrt(3.0_8)
vuc(u)=sqrt(2.0_8)
zzp(u)=zp(u)
enddo

!Selecci�n de dominio considerando los valores de calado:
!--------------------------------------------------------
 open(unit=7,file='C:\mallasubinicial.txt',status='old')
 open(2,status= 'scratch') 
 read(7,21)j
 read(7,'(A)')a
 do u=1,j
 read(7,26) no(1),no(4),no(2),no(5),no(3),no(6)
   !Selecci�n inicial. Aqu� s�lo se consideran los elementos con todos sus nodos de calado nulo.  
   if ((vt(2*i+no(1)).le.0.0).and.(vt(2*i+no(2)).le.0.0).and.(vt(2*i+no(3)).le.0.0)) then
   uu=uu+1
   write(2,27)uu,no(1),no(4),no(2),no(5),no(3),no(6)
   do uuu=1,6
   eval(no(uuu))=0.0_8
   enddo
   endif
 enddo
 rewind(7)
 read(7,23)a
 do u=1,j
 read(7,26)no(1),no(4),no(2),no(5),no(3),no(6)
   do ui=1,3
   b=ui+1-sb(ui)
   c=ui+2-sb(ui+1)      
   !Se toman los elementos que rodean a la selecci�n inicial (en ellos no se ha calculado el flujo superficial). Tambi�n se toman elementos pegados 
   !a nodos aislados o a l�neas que no forman elementos (regueros) donde no se calcul� soluci�n superficial (calados nulos). Es probable 
   !encontrar regueros pegados a la selecci�n inicial. Se trata de los elementos interfaz (definidos en la subrutina nuevamalla).
   if ((vt(2*i+no(ui)).gt.0.0_8).and.(vt(2*i+no(b)).le.0.0_8).and.(vt(2*i+no(c)).le.0.0_8)) then  
	!Elementos con un nodo apoyado en el contorno m�vil - tienen nodos esquina con: un calado positivo y dos nulos o negativos.
	uu=uu+1
    write(2,27)uu,no(1),no(4),no(2),no(5),no(3),no(6)
	do uuu=1,6
    eval(no(uuu))=0.0_8
    enddo
   elseif ((vt(2*i+no(ui)).gt.0.0_8).and.(vt(2*i+no(b)).gt.0.0_8).and.(vt(2*i+no(c)).le.0.0_8)) then
   !Elementos con dos nodos apoyados en el contorno m�vil - tienen nodos esquina con: dos calado positivos y uno nulo o negativo.
   !No hay que considerar elementos con dos calados nulos y uno positivo para dar CC pues a�n en caso de que haya dos o m�s elementos juntos de 
   !este tipo compartir�n este nodo de calado positivo, y este nodo de la malla formar� finalmente parte de un elemento de los analizados aqu�. 
	uu=uu+1
    write(2,27)uu,no(1),no(4),no(2),no(5),no(3),no(6)
	do uuu=1,6
    eval(no(uuu))=0.0_8
    enddo
	!Se llega a una interfaz sumando elementos a la selecci�n inicial. Ser� la interfaz que se calcule en la subrutina mallaaguassomeras.
	!Por tanto la interfaz ser� la misma para los dominios superficial y subterr�neo. Por tanto se pasan valores que ya tienen esos nodos (dados
	!en la subrutina nuevamalla). Ser�n CC de altura.
	!*** Se han eliminado dos l�neas que hab�a aqu�   		  	 
	!Se modifica la cota del sustrato impermeable para garantizar que el flujo que entra o sale por el contorno m�vil para el modelo subterr�neo
	!sea el que sale o entra por el contorno m�vil para el modelo superficial. 
	!*** zzp s�lo se usa para seleccionar el contorno m�vil
	zzp(no(ui))=z(no(ui))
	zzp(no(b))=z(no(b))	      
   endif
   enddo   
 enddo

!Arreglo de las condiciones de contorno:
!---------------------------------------
!Fuera de la nueva malla considerada (en todo el dominio que no es subterr�neo) se dan condiciones de contorno nulas 
!de nivel fre�tico (hd). As� no se calcular� la soluci�n en estos nodos.	
do u=1,i
 if (eval(u).eq.sqrt(3.0_8)) then	   
 vuc(u)=0.0_8
 endif
enddo
	
!Introducci�n de las CC de flujo subterr�neo que existen en el contorno no m�vil (le�das desde fichero).    
do u=1,i+2
 read(7,'(A)')a
enddo
eaf=0
do while (eaf.eq.0)
 read (7,46,iostat=eaf) u,fe,fi
 io=io+1
 
!En la intersecci�n del contorno m�vil con el contorno del dominio conjunto tambi�n se tienen las CC de flujo subterr�neo para el contorno del 
!dominio conjunto. Los nodos de dicha intersecci�n pertenecer�n a la zona de CC de flujo superficial del dominio conjunto (el contorno m�vil est� 
!en la orilla del dominio superficial). As�, en el fichero para flujo subterr�neas se tiene CC de velocidad nula en ellos. 
!Por este motivo s�lo se utiliza el valor dado de altura de la l�mina en la intersecci�n.
if (eaf.eq.0) then  
  if (io.le.i) then
	if ((vuc(u).eq.sqrt(2.0_8)).and.(fe.ne.'o')) then    
	qx(u)=fi
	endif
  elseif ((io.gt.i).and.(io.le.2*i)) then	 
	if ((vuc(u).eq.sqrt(2.0_8)).and.(fe.ne.'o')) then    
	qy(u)=fi
	endif															   
  else
    !En el contorno m�vil, vuc ser� distinto de sqrt(2.0) ya que se impuso CC. As�, no se seleccionan los nodos de la intersecci�n.
	if ((vuc(u).eq.sqrt(2.0_8)).and.(fe.ne.'o')) then 	   
	vuc(u)=fi
	endif
  endif
 endif 
enddo
		
!Escritura de la malla:
!----------------------
!Se escribe una nueva malla para subterr�neo que no se va a modificar durante la resoluci�n de la ecuaci�n de agua subterr�nea.
!La malla tendr� la i original porque los n�meros de nodo ir�n hasta i, y la nueva j para construir bien las cajas existentes.
rewind(2)                             
close(7)           
open(3,file = 'C:\mallasub.txt',status='replace')
write(3,'(A)')'nudos'
write(3,10)i
write(3,'(A)')'elementos'
write(3,10)uu
write(3,'(A)')'elemento p1  p2  p3  p4  p5  p6 (antihorarios)'
do u=1,uu
 read(2,'(A)')a
 write(3,'(A)')a
enddo
close(2)
close(3)
j=uu
end

!----------------------------------------------------------------------------------------------------------------------------------
!Subrutinas para el c�lculo de las dimensiones de los vectores que llevan la matriz del sistema (predimensionamiento).
!----------------------------------------------------------------------------------------------------------------------------------
!Subrutina DIMVECTAS
!C�lculo de la dimensi�n para el sistema de ecuaciones para flujo superficial.
!Se consideran 9 posiciones en la matriz del sistema. La primera posici�n ir�a desde la fila 1 a la 'i' y de la columna 1 a la 'i'. 
!----------------------------------------------------------------------------------------------------------------------------------		  
subroutine dimvectas (i,j,newton,sino,vn,est)
use elemental
use allocatacion
integer*4, dimension(:),allocatable::nt,Eas,Easdosi,Eastresi,Eeqe,Eqqe,Eqe,Eee
integer*4 i,j,vn(i),Et
character newton*2,sino*2,est*2 

allocate(nt(i),Eas(3*i),Easdosi(3*i),Eastresi(3*i),Eeqe(i),Eqqe(i),Eqe(i),Eee(i))

39     format(4/,A80)
40     format(6X,6(X,I5))

!Inizializaci�n previa de variables:
!-----------------------------------
Et=0
do u=1,i
nt(u)=0
Eqqe(u)=0
Eqe(u)=0
Eeqe(u)=0
Eee(u)=0
enddo

!C�lculo de la dimensi�n necesaria sin ensamblar para una caja determinada
!-------------------------------------------------------------------------
open(unit=1,file='C:\malla.txt',status='old')
read(1,39)ac
do u=1,j
read(1,40) n(1),n(4),n(2),n(5),n(3),n(6) 
 !C�lculo en 'nt(i)' del n�mero de elementos que tiene el nodo 'i' en com�n.
 do ui=1,6
 nt(n(ui))=nt(n(ui))+1
 enddo
enddo
close(1)

!Se calcula el n�mero de coeficientes por fila que se podr�an generar. 
!'Eeqe' ser� el n� de coeficientes en los nodos esquina si se eval�a una funci�n en nodos esquina y cuadr�ticos, 'Eqqe' el n� en los 
!nodos cuadr�ticos si se eval�a en nodos esquina y cuadr�ticos, 'Eee' el n� en los nodos esquina si se eval�a en los nodos esquina.
!'Eqe' el n� en los nodos cuadr�ticos si se eval�a en los nodos esquina.   
do u=1,i
 if (vn(u).eq.0) then	
 !Si 'u' es nodo esquina. 'Eqe(u)' y 'Eqqe(u)' ser�n nulos. 
  Eeqe(u)=6*nt(u)
  Eee(u)=3*nt(u) 
 else
 !Si 'u' no es nodo esquina. 'Eee(u)' y 'Eeqe(u)' ser�n nulos.					 
  Eqqe(u)=6*nt(u)
  Eqe(u)=3*nt(u)
 endif
enddo
!As� se pueden conocer el n� de coeficientes por fila para una caja determinada. S�lo se usar�a Eeqe si se tuviese 
!una caja con funciones peso discretizadas en elemento lineal y funci�n a integrar en elemento cuadr�tico.

!C�lculo de la dimensi�n para la matriz del sistema suponiendo que s�lo se a�ada una caja all� donde haya cajas no nulas 
!-----------------------------------------------------------------------------------------------------------------------
!S�lo ser�a necesario calcular la variable 'Eas'.
if ((newton.eq.'no').and.(est.eq.'no')) then
!En este caso en las posiciones de la matriz del sistema (1,2) y (2,1) no se calcula ninguna caja. 
 do u=1,i
 !Se tienen en cuenta los tipos de cajas que hay y se calcula el n� total de coeficientes por fila.
 !En 'Eas(u)' est� el n� de coeficientes por fila 'u'. Desde la fila 1 a la 'i':
 Eas(u)=Eeqe(u)+Eqqe(u)+Eee(u)+Eqe(u)
 !Desde la fila 'i+1' a la '2*i':
 Eas(u+i)=Eeqe(u)+Eqqe(u)+Eee(u)+Eqe(u)
  !Desde la fila '2*i+1' a la '3*i': 
  if (sino.eq.'si')then
  Eas(u+2*i)=2*Eeqe(u)+Eee(u)
  else
  !No se calcula caja para la posici�n (3,3).
  Eas(u+2*i)=2*Eeqe(u)
  endif
 !Se tienen en cuenta los tipos de cajas que hay y se calcula el n� de coeficientes hasta la columna 'i' por fila.
 Easdosi(u)=Eeqe(u)+Eqqe(u)
 Easdosi(u+i)=0
 Easdosi(u+2*i)=Eeqe(u)
 !Se tienen en cuenta los tipos de cajas que hay y se calcula el n� de coeficientes hasta la columna '2*i' por fila.
 Eastresi(u)=Eeqe(u)+Eqqe(u)
 Eastresi(u+i)=Eeqe(u)+Eqqe(u)
 Eastresi(u+2*i)=2*Eeqe(u)
 enddo
else	
!En otro caso se considera que hay cajas en todas las posiciones.
 do u=1,i
 !Se calcula el n� total de coeficientes por fila.
 Eas(u)=2*Eeqe(u)+2*Eqqe(u)+Eee(u)+Eqe(u)
 Eas(u+i)=2*Eeqe(u)+2*Eqqe(u)+Eee(u)+Eqe(u)
 Eas(u+2*i)=2*Eeqe(u)+Eee(u)
 !Se calcula el n� de coeficientes hasta la columna 'i' por fila.
 Easdosi(u)=Eeqe(u)+Eqqe(u)
 Easdosi(u+i)=Eeqe(u)+Eqqe(u)
 Easdosi(u+2*i)=Eeqe(u)
 !Se calcula el n� de coeficientes hasta la columna '2*i' por fila.
 Eastresi(u)=2*Eeqe(u)+2*Eqqe(u)
 Eastresi(u+i)=2*Eeqe(u)+2*Eqqe(u)
 Eastresi(u+2*i)=2*Eeqe(u)
 enddo
endif
!C�lculo de la dimensi�n total que tendr�n los vectores. Es necesario reservar un espacio de '3*i+1' para referenciar 
!en que posici�n empiezan los coeficientes de cada fila (se hace a continuaci�n).
do u=1,i
Et=Et+Eas(u)+Eas(u+i)+Eas(u+2*i)
enddo
Et=Et+3*i+1

!Se referencia el vector puntero
!-------------------------------
allocate (cia(Et),ca(Et),ita(Et),sa(Et))
do u=1,Et
ita(u)=0
sa(u)=0.0_8
enddo
!Se referencia en cada componente 'u' del vector en que posici�n se empezar�n a escribir los coeficientes de la fila 'u'.
!La suma del n� total de coeficientes por fila 'Eas' dar� informaci�n de la posici�n del primer coeficiente de la siguiente fila.  
!Finalmente 'ita(3*i+1)'-1 es igual a 'Et' (posteriormente ndim).
ita(1)=3*i+2	   
do u=1,3*i
ita(u+1)=ita(u)+Eas(u)
enddo

!C�lculo de otras variables para no reservar m�s coeficientes para cajas que se sumen sobre otras.
!-------------------------------------------------------------------------------------------------
!Se calcula la variable 'posi(u)' con posici�n ocupada en 'ita' por el primer coeficiente de cada fila 'u'. Tambi�n la variable 
!'posdosi(u)' con posici�n ocupada en 'ita' por el primer coeficiente situado en una columna posterior a 'i' de cada fila 'u'. 
!La suma del n� de coeficientes por fila 'Easdosi' dar� esta informaci�n. Tambi�n la variable 'postresi(u)' con posici�n ocupada 
!en 'ita' por el primer coeficiente situado en una columna posterior a '2*i' de cada fila 'u'.
do u=1,3*i
posi(u)=ita(u)
posdosi(u)=ita(u)+Easdosi(u)
postresi(u)=ita(u)+Eastresi(u)
enddo
!Al escribir cajas en mismos lugares basta con inicializar posiciones utilizando estas variables y sumar los nuevos coeficientes 
!a los existentes. Esto se puede hacer dado que se crean coeficientes cuya referencia sigue el mismo orden (los elementos se leen 
!por el mismo orden). 

deallocate(nt,Eas,Easdosi,Eastresi,Eeqe,Eqqe,Eqe,Eee)
end

!---------------------------------------------------------------------------------------------------------------------------------
!Subrutina DIMVECTSB
!C�lculo de la dimensi�n para el sistema de ecuaciones para flujo subterr�neo.
!Se considera que solo hay 1 posici�n en la matriz del sistema. 
!---------------------------------------------------------------------------------------------------------------------------------
subroutine dimvectsb (i,j,vn)
use elemental
use allocatacion
integer*4, dimension(:),allocatable::nt,Esb,Eee
integer*4 i,j,vn(i),Et	

allocate(nt(i),Esb(i),Eee(i))

39     format(4/,A80)
40     format(6X,6(X,I5))

!Inizializaci�n previa de variables:
!-----------------------------------
Et=0
do u=1,i
nt(u)=0
Eee(u)=0
enddo

!C�lculo de la dimensi�n necesaria sin ensamblar para una caja determinada
!-------------------------------------------------------------------------
open(unit=1,file='C:\mallasub.txt',status='old')
read(1,39)ac
do u=1,j
read(1,40) n(1),n(4),n(2),n(5),n(3),n(6) 
 do ui=1,6
 nt(n(ui))=nt(n(ui))+1
 enddo
enddo
close(1)
    
do u=1,i
 if (vn(u).eq.0) then
 !Si 'u' es nodo esquina. Eee(u) ser� nulo si 'u' no es esquina.	 
  Eee(u)=3*nt(u) 
 endif
enddo
!S�lo habr� caja con funciones peso discretizadas en elemento lineal y funci�n a integrar en elemento lineal.

!C�lculo de la dimensi�n para la matriz del sistema suponiendo que s�lo se a�ada una caja 
!----------------------------------------------------------------------------------------
do u=1,i
Esb(u)=Eee(u)
enddo
do u=1,i
Et=Et+Esb(u)
enddo
Et=Et+i+1

!Se referencia el vector puntero
!-------------------------------
allocate (cia(Et),ca(Et),ita(Et),sa(Et))
do u=1,Et
ita(u)=0
sa(u)=0.0_8
enddo
ita(1)=i+2	   
do u=1,i
ita(u+1)=ita(u)+Esb(u)
enddo

!C�lculo de otras variables para no reservar m�s coeficientes para cajas que se sumen sobre otras.
!-------------------------------------------------------------------------------------------------
!Se calcula una variable 'pos(u)' con posici�n ocupada en 'ita' por el primer coeficiente de cada fila 'u'.
do u=1,i
pos(u)=ita(u)
enddo

deallocate(nt,Esb,Eee)
end

!----------------------------------------------------------------------------------------------------------------------------------
!Subrutinas para el c�lculo de coeficientes de las integrales de contorno (vectores elementales).
!----------------------------------------------------------------------------------------------------------------------------------
!Subrutina VECTORCONTORNOPRESIONES
!En esta subrutina se calculan las integrales de contorno para la ecuaci�n de Navier-Stokes 2D o para la ecuaci�n de aguas someras.
!La integraci�n se ha hecho anal�ticamente (a mano) obteniendo las expresiones necesarias. 
!----------------------------------------------------------------------------------------------------------------------------------
subroutine vectorcontornopresiones(i,j,x,y,vv,vvv,vectorx,vectory,vector,nu,ten)
integer*4 i,j,u,ui,sb(4),n(6),a,b,c,d,e,f
real*8 vectorx(i),vectory(i),vector(i),vv(3*i),vvv(3*i),x(i),y(i),xa,xb,ya,yb,pu,pv,qu,qv,nu,ax,ay,ma,mb,mc,md,jac            
character ac*80,ten*2
 
39     format(4/,A80)
40     format(6X,6(X,I5))

!Inizializaci�n previa de variables:
!-----------------------------------
ma=3.0_8/20.0_8
mb=1.0_8/60.0_8
mc=1.0_8/5.0_8
md=2.0_8/15.0_8
sb=(/0,0,3,3/)

!C�lculo de los valores para cada lado del elemento
!-------------------------------------------------- 
!Ensamblaje directo al sobreescribir los vectores con la suma de su valor y el de los coeficientes generados.
open(unit=1,file='C:\malla.txt',status='old')
read(1,39)ac
do u=1,j
!Selecci�n de cada elemento
read(1,40) n(1),n(4),n(2),n(5),n(3),n(6)
jac=abs((x(n(1))-x(n(3)))*(y(n(2))-y(n(3)))-(x(n(2))-x(n(3)))*(y(n(1))-y(n(3))))
 do ui=1,3
  !Selecci�n de cada lado del elemento.
  a=ui
  b=ui+1-sb(ui)
  c=ui+2-sb(ui+1)
  d=ui+3
  e=ui+4-sb(ui)
  f=ui+5-sb(ui+1)
 
 !C�lculo previo de coeficientes necesarios
 !-----------------------------------------
 !Coeficientes calculados por separado por la longitud de la expresi�n.     	
 pu=(vv(n(c))*ma-vv(n(b))*mb+vv(n(e))*mc)*vv(2*i+n(c))+(vv(n(c))*mb+vv(n(b))*mb+vv(n(e))*md)*vv(2*i+n(b))
 pv=(vv(i+n(c))*ma-vv(i+n(b))*mb+vv(i+n(e))*mc)*vv(2*i+n(c))+(vv(i+n(c))*mb+vv(i+n(b))*mb+vv(i+n(e))*md)*vv(2*i+n(b)) 
 qu=(vv(n(c))*mb+vv(n(b))*mb+vv(n(e))*md)*vv(2*i+n(c))+(-vv(n(c))*mb+vv(n(b))*ma+vv(n(e))*mc)*vv(2*i+n(b))
 qv=(vv(i+n(c))*mb+vv(i+n(b))*mb+vv(i+n(e))*md)*vv(2*i+n(c))+(-vv(i+n(c))*mb+vv(i+n(b))*ma+vv(i+n(e))*mc)*vv(2*i+n(b))
 !C�lculo de los catetos correspondientes al producto de la longitud del lado por coseno o seno del �ngulo que forman la direcci�n normal 
 !y las direcciones de los ejes. 
 xa=x(n(c))-x(n(b)) 
 xb=x(n(a))-x(n(b))		   
 ya=y(n(c))-y(n(b))		
 yb=y(n(a))-y(n(b))
 !C�lculo del signo de las componentes del vector normal al contorno que tiene orientaci�n saliente.
 call direccion(x,y,i,n(b),n(c),n(a),ax,ay)
 
 !Para integrales de contorno de las ecuaciones din�micas - t�rmino de presi�n
 !----------------------------------------------------------------------------
 vectorx(n(c))=vectorx(n(c))-(vvv(2*i+n(c))/6.0_8)*ax*abs(ya)*9.81_8 
 vectorx(n(b))=vectorx(n(b))-(vvv(2*i+n(b))/6.0_8)*ax*abs(ya)*9.81_8
 vectorx(n(e))=vectorx(n(e))-((vvv(2*i+n(b))+vvv(2*i+n(c)))/3.0_8)*ax*abs(ya)*9.81_8 
 vectory(n(c))=vectory(n(c))-(vvv(2*i+n(c))/6.0_8)*ay*abs(xa)*9.81_8 
 vectory(n(b))=vectory(n(b))-(vvv(2*i+n(b))/6.0_8)*ay*abs(xa)*9.81_8 
 vectory(n(e))=vectory(n(e))-((vvv(2*i+n(b))+vvv(2*i+n(c)))/3.0_8)*ay*abs(xa)*9.81_8 
 
 !Para integrales de contorno de las ecuaciones din�micas - t�rmino viscoso	(tensiones)
 !-------------------------------------------------------------------------------------
 !El usuario decide si considerar estos integrales o no. Tras ensamblar generar�n valores no nulos en los nodos internos de la malla
 !dado que no existe continuidad de la derivada de las funciones de interpolaci�n en la direcci�n normal (s�lo la habr�a en la direcci�n del lado).
 !As�, estos t�rminos pueden dar lugar a una mayor inestabilidad.
 if (ten.eq.'si') then
 vectorx(n(c))=vectorx(n(c))+(3.0_8*vv(n(c))+vv(n(b))-4.0_8*vv(n(e)))*(yb*ax*abs(ya)-xb*ay*abs(xa))*nu/(6.0_8*jac)
 vectorx(n(b))=vectorx(n(b))+(-vv(n(c))-3.0_8*vv(n(b))+4.0_8*vv(n(e)))*(yb*ax*abs(ya)-xb*ay*abs(xa))*nu/(6.0_8*jac)
 vectorx(n(e))=vectorx(n(e))+4.0_8*(vv(n(c))-vv(n(b)))*(yb*ax*abs(ya)-xb*ay*abs(xa))*nu/(6.0_8*jac)
 vectory(n(c))=vectory(n(c))+(3.0_8*vv(i+n(c))+vv(i+n(b))-4.0_8*vv(i+n(e)))*(yb*ax*abs(ya)-xb*ay*abs(xa))*nu/(6.0_8*jac)
 vectory(n(b))=vectory(n(b))+(-vv(i+n(c))-3.0_8*vv(i+n(b))+4.0_8*vv(i+n(e)))*(yb*ax*abs(ya)-xb*ay*abs(xa))*nu/(6.0_8*jac)
 vectory(n(e))=vectory(n(e))+4.0_8*(vv(i+n(c))-vv(i+n(b)))*(yb*ax*abs(ya)-xb*ay*abs(xa))*nu/(6.0_8*jac)
 endif

 !Para integrales de contorno de la ecuaci�n de continuidad
 !---------------------------------------------------------
 !pu=q referente al nodo n(c) de integrar el calado por u. pv=q referente al nodo n(c) de integrar el calado por v
 !qu=q referente al nodo n(b) de integrar el calado por u. qv=q referente al nodo n(b) de integrar el calado por v
 vector(n(c))=vector(n(c))+(pu*ax*abs(ya)+pv*ay*abs(xa))
 vector(n(b))=vector(n(b))+(qu*ax*abs(ya)+qv*ay*abs(xa))  
	
 enddo
enddo
close(1)
end

!--------------------------------------------------------------------------------------------------
!Subrutina VECTORCONTORNOCAUDALES
!En esta subrutina se calculan las integrales de contorno para la ecuaci�n para flujo subterraneo.
!La integraci�n se ha hecho anal�ticamente (a mano) obteniendo las expresiones necesarias.
!--------------------------------------------------------------------------------------------------
subroutine vectorcontornocaudales(i,j,x,y,zp,zzp,vector,qxx,qyy,vb,al)         !*** A�adidos zp,zzp,vb,al
use elemental								                                   !*** A�adida esta l�nea
integer*4, dimension(:),allocatable::po                                        !*** A�adida esta l�nea
integer*4 i,j,sb(4),a,b,c										               !*** Eliminados n(6),u,ui
real*8 vector(i),x(i),y(i),ax,ay,qxx(i),qyy(i),zp(i),zzp(i),al,vb(3*i),Nn(3,3) !*** Eliminados xa,ya y a�adidos zp(i),zzp(i),al,vb(3*i),Nn(3,3)
																               !*** Eliminado character ac*80

allocate(po(i))    !*** A�adida esta l�nea 

39     format(4/,A80)
40     format(6X,6(X,I5))
		  
!Inizializaci�n previa de variables:
!-----------------------------------
sb=(/0,0,3,3/)
do u=1,i		   !*** A�adida esta l�nea
po(u)=pos(u)	   !*** A�adida esta l�nea
enddo			   !*** A�adida esta l�nea

!C�lculo de los valores para cada lado del elemento
!--------------------------------------------------
open(unit=3,file='C:\mallasub.txt',status='old')
read(3,39)ac
do u=1,j
read(3,40)n(1),n(4),n(2),n(5),n(3),n(6) 
 do uu=1,3		   !*** A�adida esta l�nea
 do uuu=1,3		   !*** A�adida esta l�nea
 Nn(uu,uuu)=0.0_8  !*** A�adida esta l�nea
 enddo			   !*** A�adida esta l�nea
 enddo			   !*** A�adida esta l�nea
 do ui=1,3
   a=ui
   b=ui+1-sb(ui)
   c=ui+2-sb(ui+1)
   xa=x(n(c))-x(n(b)) 
   ya=y(n(c))-y(n(b)) 
   call direccion(x,y,i,n(b),n(c),n(a),ax,ay)
   
   if ((zp(n(c)).eq.zzp(n(c))).or.(zp(n(b)).eq.zzp(n(b)))) then !*** A�adida esta l�nea
   !C�lculo de las Integrales de contorno de caudal	(valor interpolado linealmente)	que aparecen en la ecuaci�n. 
   vector(n(c))=vector(n(c))+((2.0_8*qxx(n(c))+qxx(n(b)))*ax*abs(ya)+(2.0_8*qyy(n(c))+qyy(n(b)))*ay*abs(xa))/6.0_8
   vector(n(b))=vector(n(b))+((2.0_8*qxx(n(b))+qxx(n(c)))*ax*abs(ya)+(2.0_8*qyy(n(b))+qyy(n(c)))*ay*abs(xa))/6.0_8
   else																																 !*** A�adida esta l�nea
   !*** Integrales de contorno diferentes para condici�n de goteo																		 
   !*** vector(n(c))=vector(n(c))-al*((2.0_8*vb(2*i+n(c))+vb(2*i+n(b)))*ax*abs(ya)+(2.0_8*vb(2*i+n(c))+vb(2*i+n(b)))*ay*abs(xa))/6.0_8	 
   !*** vector(n(b))=vector(n(b))-al*((2.0_8*vb(2*i+n(b))+vb(2*i+n(c)))*ax*abs(ya)+(2.0_8*vb(2*i+n(b))+vb(2*i+n(c)))*ay*abs(xa))/6.0_8   
   !*** Nn(c,c)=-al*2.0_8*(ax*abs(ya)+ay*abs(xa))/6.0_8																					 
   !*** Nn(c,b)=-al*1.0_8*(ax*abs(ya)+ay*abs(xa))/6.0_8																					 
   !*** Nn(b,c)=-al*1.0_8*(ax*abs(ya)+ay*abs(xa))/6.0_8																					 
   !*** Nn(b,b)=-al*2.0_8*(ax*abs(ya)+ay*abs(xa))/6.0_8																					 
   vector(n(c))=vector(n(c))-al*(2.0_8*vb(2*i+n(c))+vb(2*i+n(b)))/6.0_8																 !*** A�adida esta l�nea
   vector(n(b))=vector(n(b))-al*(2.0_8*vb(2*i+n(b))+vb(2*i+n(c)))/6.0_8  															 !*** A�adida esta l�nea
   Nn(c,c)=al*2.0_8/6.0_8																											 !*** A�adida esta l�nea
   Nn(c,b)=al*1.0_8/6.0_8																											 !*** A�adida esta l�nea
   Nn(b,c)=al*1.0_8/6.0_8																											 !*** A�adida esta l�nea
   Nn(b,b)=al*2.0_8/6.0_8																											 !*** A�adida esta l�nea
   endif																															 !*** A�adida esta l�nea
   !*** Cambio el comentario de dos l�neas que hab�a por el de debajo
   !*** En el contorno m�vil las integrales de contorno se calculan con caudales que no dependen de velx y vely sino de la diferencia 
   !*** de altura entre flujo superficial y subterr�neo (la velocidad es impl�cita en la f�rmula). Estos valores calculados se tienen 
   !*** en cuenta al no imponer condiciones de nivel sino de caudal en el contorno m�vil.
 enddo
 !*** Necesario hacer sumas para cada elemento para rellenar todos los espacios en su posici�n.
 !*** S�lo habr� un lado que genera coeficientes (si hubiese m�s el elemento pertenece a As) y por tanto no hace falta hacer sumas en Nn. Se 
 !*** podr�a poner un exit arriba.
 call suma(Nn,n,0,0,i,3,3,po)	!*** A�adida esta l�nea
 enddo
close(3)
deallocate(po)	!*** A�adida esta l�nea
end

!------------------------------------------------------------------------------------------------------------------------------------
!Subrutina DIRECCION
!Esta subrutina obliga a que la orientaci�n del vector normal en cada lado del contorno elemental sea saliente respecto al elemento.
!Cada vez que se utiliza se analiza un lado del elemento. As�, se analiza el lado formado por los nodos 'k' y 'm'.
!'a' es el signo de la componente horizontal del vector, 'b' es el signo de la componente vertical del vector.
!------------------------------------------------------------------------------------------------------------------------------------
subroutine direccion(x,y,i,k,m,p,a,b)
integer*4 i,k,m,p
real*8 x(i),y(i),a,b

!La direcci�n del lado es distinta a la horizontal o a la vertical.
if (((x(k)-x(m)).ne.0.0_8).and.((y(k)-y(m)).ne.0.0_8)) then
  if ((y(p)-y(k)-((y(k)-y(m))/(x(k)-x(m)))*(x(p)-x(k))).lt.0.0_8) then
  b=1.0_8
	if ((x(p)-x(k)-((x(k)-x(m))/(y(k)-y(m)))*(y(p)-y(k))).lt.0.0_8) then
	a=1.0_8
	elseif ((x(p)-x(k)-((x(k)-x(m))/(y(k)-y(m)))*(y(p)-y(k))).gt.0.0_8) then
	a=-1.0_8
	endif
  elseif ((y(p)-y(k)-((y(k)-y(m))/(x(k)-x(m)))*(x(p)-x(k))).gt.0.0_8) then
  b=-1.0_8
    if ((x(p)-x(k)-((x(k)-x(m))/(y(k)-y(m)))*(y(p)-y(k))).lt.0.0_8) then
    a=1.0_8
	elseif ((x(p)-x(k)-((x(k)-x(m))/(y(k)-y(m)))*(y(p)-y(k))).gt.0.0_8) then
	a=-1.0_8
	endif
  endif
!La direcci�n del lado es vertical.
elseif ((x(k)-x(m)).eq.0.0_8) then
b=0.0_8 
  if ((x(k)-x(p)).gt.0.0_8) then
  a=1.0_8
  elseif ((x(k)-x(p)).lt.0.0_8) then
  a=-1.0_8
  endif 
!La direcci�n del lado es horizontal.
elseif ((y(k)-y(m)).eq.0.0_8) then
a=0.0_8 
  if ((y(p)-y(k)).lt.0.0_8) then
  b=1.0_8
  elseif ((y(p)-y(k)).gt.0.0_8) then
  b=-1.0_8
  endif
endif 				 
end

!-----------------------------------------------------------------------------------------------------------------------------------------
!Subrutinas para el c�lculo de coeficientes de las integrales en el dominio (formulaci�n Bubnov-Galerkin).
!Las funciones de forma est�n escrita en coordenadas locales. 
!Para elementos lineales 'Mp' ser�n las funciones de forma, 'Mpx' las derivadas en x de las funciones de forma, 'Mpy' las derivadas en y 
!de las funciones de forma. Para elementos cuadr�ticos 'Mi' ser�n las funciones de forma, 'Mix' ser�n las derivadas en x de las funciones 
!de forma, 'Miy' las derivadas en y de las funciones de forma.
!-----------------------------------------------------------------------------------------------------------------------------------------
!Subrutina VECTORCONTORNOFUENTE
!En esta subrutina se imponen las condiciones interiores de lluvia para flujo superficial calculando vectores elementales 
!que ir�n al t�rmino independiente del sistema.	
!Se hace una integraci�n en el dominio de la intensidad de lluvia (tiene unidades m3/m2*s), cuyos valores nodales est�n guardados en 
!la variable 'ql'. Las integrales son calculadas num�ricamente mediante cuadraturas.
!La funci�n intensidad se interpola s�lo a trav�s de los nodos esquina (interpolaci�n lineal). Es irrelevante como se interpole ya que 
!los valores nodales ser�n id�nticos e iguales a la intensidad que define el usuario para todo el dominio. 
!'n' es una variable global, tiene dimensi�n 6 y lleva los n�meros de los nodos de cada elemento le�dos desde fichero. 
!-----------------------------------------------------------------------------------------------------------------------------------------
subroutine vectorcontornofuente(i,j,x,y,vectorz,ql)
use elemental
integer*4 i,j
real*8 x(i),y(i),Ip,Mp(3),vectorz(i),ql(i)            

39     format(4/,A80)
40     format(6X,6(X,I5))

!Inizializaci�n previa de variables:
!-----------------------------------
!En 'luno,ldos,ja', que son variables globales, van los puntos y pesos de integraci�n para la cuadratura.
!Ciertos valores son nulos dado que su dimensi�n es 7 y s�lo se usar�n 3 puntos de integraci�n. 
luno=(/0.5_8, 0.0_8, 0.5_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8/)		   
ldos=(/0.5_8, 0.5_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8/)
ja=(/1.0_8/6.0_8, 1.0_8/6.0_8, 1.0_8/6.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8/)

!C�lculo de los valores para cada elemento
!-----------------------------------------
!Se calculan valores s�lo para elementos del dominio, que est�n guardados en malla.txt
open(unit=1,file='C:\malla.txt',status='old')
read(1,39)ac
 do u=1,j
 !Los nodos en el fichero est�n escritos teniendo en cuenta el sentido antihorario. Las funciones de forma para elementos lineales 
 !est�n definidas de forma que est�n referidas a los nodos esquina en sentido antihorario. Por ello se almacenan los nodos en este orden.
 read(1,40) n(1),n(4),n(2),n(5),n(3),n(6)

	 !Se calculan antes estos valores ya que son contantes para cada elemento. Siempre se calcular�a lo mismo si se eval�an estas expresiones 
	 !para cada punto de integraci�n. Dar�a igual meterlo en el siguiente bucle aunque conllevar�a mayor coste computacional.
	 xa=x(n(1))-x(n(3)) 
	 xb=x(n(2))-x(n(3)) 
	 ya=y(n(1))-y(n(3))
	 yb=y(n(2))-y(n(3))  
	 jac=abs(xa*yb-xb*ya)	

	 !Integraci�n con tres puntos. Polinomios de grado superior a grado 2 quedar�an integrados con error.
	 do uu=1,3	      
	 !Funciones necesarias.
	 Mp(1)=luno(uu)
	 Mp(2)=ldos(uu)
	 Mp(3)=1.0_8-luno(uu)-ldos(uu)	 
	 
	  !Se conocen valores de intensidad en los nodos esquina.
	  Ip=Mp(1)*ql(n(1))+Mp(2)*ql(n(2))+Mp(3)*ql(n(3))
	  !Ensamblaje directo de los vectores elementales.
	  do uuu=1,3
	  vectorz(n(uuu))=vectorz(n(uuu))+(Mp(uuu)*Ip)*ja(uu)*jac 	  
	  enddo	
	 enddo
 enddo
close(1)
end

!---------------------------------------------------------------------------------------------------------------------------------------
!Subrutina VECTORCONTORNOFUENTESUB
!En esta subrutina se imponen las condiciones interiores de lluvia para flujo subterr�neo calculando vectores elementales 
!que ir�n al t�rmino independiente del sistema.	
!Se hace una integraci�n en el dominio de la intensidad de lluvia (tiene unidades m3/m2*s), cuyos valores nodales est�n guardados en 
!la variable 'ql'. 
!La funci�n intensidad se interpola s�lo a trav�s de los nodos esquina (interpolaci�n lineal). 
!---------------------------------------------------------------------------------------------------------------------------------------
subroutine vectorcontornofuentesub(i,j,x,y,vectorz,ql)
use elemental
integer*4 i,j
real*8 x(i),y(i),Ip,Mp(3),vectorz(i),ql(i)            

39     format(4/,A80)
40     format(6X,6(X,I5))

!Inizializaci�n previa de variables:
!-----------------------------------
luno=(/0.5_8, 0.0_8, 0.5_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8/)		   
ldos=(/0.5_8, 0.5_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8/)
ja=(/1.0_8/6.0_8, 1.0_8/6.0_8, 1.0_8/6.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8/)

!C�lculo de los valores para cada elemento
!-----------------------------------------
!Se calculan valores s�lo para elementos del dominio, que est�n guardados en mallasub.txt
open(unit=3,file='C:\mallasub.txt',status='old')
read(3,39)ac
 do u=1,j
 read(3,40) n(1),n(4),n(2),n(5),n(3),n(6)

	 xa=x(n(1))-x(n(3)) 
	 xb=x(n(2))-x(n(3)) 
	 ya=y(n(1))-y(n(3))
	 yb=y(n(2))-y(n(3))  
	 jac=abs(xa*yb-xb*ya)	

	 do uu=1,3		      
	 
	 Mp(1)=luno(uu)
	 Mp(2)=ldos(uu)
	 Mp(3)=1.0_8-luno(uu)-ldos(uu)
	 
	  !Se conocen valores de intensidad en los nodos esquina.
	  Ip=Mp(1)*ql(n(1))+Mp(2)*ql(n(2))+Mp(3)*ql(n(3))
	  do uuu=1,3
	  vectorz(n(uuu))=vectorz(n(uuu))+(Mp(uuu)*Ip)*ja(uu)*jac 	  
	  enddo	
	 enddo
 enddo
close(3)
end

!-----------------------------------------------------------------------------------------------------------------------------------------
!Subrutina CAJASAB
!A continuaci�n se calculan coeficientes de las matrices elementales que forman las caja A, Bx y By.
!'Nn' son las matrices elementales de la caja A, 'Nm' son las matrices elementales de la caja Bx y 'Nk' son las matrices elementales 
!de la caja By. Aparecen en las ecuaciones din�micas para flujo superficial.
!'sa' es el vector que lleva los coeficientes de la matriz del sistema no nulos de fuera de la diagonal y 'ita' es el vector puntero 
!que lleva su posici�n (a�n sin formato MSR). Son variables globales que ya han sido definidas.
!'poi,podosi,potresi' son inicializadas con las variables globales 'posi,posdosi,postresi' para indicar en que posici�n del vector para 
!cada fila se empiezan a escribir los coeficientes. De esta forma, al escribir los coeficientes uno tras otro, se suman de forma 
!adecuada si se escriben dos cajas de matrices en la misma posici�n de la matriz del sistema.
!Por ejemplo 'poi' lo indica para las posiciones (1,1),(2,1) y (3,1); 'podosi' lo indica para las posiciones (1,2),(2,2) y (3,2);... 
!Los coeficientes de las matrices elementales se integran y se escriben sobre los vectores ita y sa a trav�s de la subrutina suma.
!----------------------------------------------------------------------------------------------------------------------------------------- 
subroutine cajasab(i,j,x,y,nu,del)
use elemental
integer*4, dimension(:),allocatable::poi,podosi,potresi
integer*4 i,j
real*8 x(i),y(i),Nn(6,6),Nm(6,3),Nk(6,3),Mix(6),Miy(6),Mp(3),nu,del

allocate(poi(3*i),podosi(3*i),potresi(3*i))
   
39     format(4/,A80)
40     format(6X,6(X,I5))

!Inizializaci�n previa de variables:
!-----------------------------------
luno=(/0.5_8, 0.0_8, 0.5_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8/)		   
ldos=(/0.5_8, 0.5_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8/)
ja=(/1.0_8/6.0_8, 1.0_8/6.0_8, 1.0_8/6.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8/)
do u=1,3*i
poi(u)=posi(u)
podosi(u)=posdosi(u)
potresi(u)=postresi(u)
enddo

!C�lculo de los valores para cada elemento
!-----------------------------------------
open(unit=1,file='C:\malla.txt',status='old')
read(1,39)ac
 do u=1,j
     !Los nodos en el fichero est�n escritos teniendo en cuenta el sentido antihorario.
	 !Las funciones de forma para elementos cuadr�ticos est�n definidas de forma que est�n referidas primero a los nodos esquina en sentido 
	 !antihorario y despu�s a los tres nodos cuadr�ticos en sentido antihorario. Por ello se almacenan los nodos en este orden.	 
     read(1,40) n(1),n(4),n(2),n(5),n(3),n(6)
  
     xa=x(n(1))-x(n(3)) 
	 xb=x(n(2))-x(n(3)) 
	 ya=y(n(1))-y(n(3))
	 yb=y(n(2))-y(n(3))
	 jac=abs(xa*yb-xb*ya)

	 !Inicializaci�n de las matrices elementales para el c�lculo de nuevas matrices elementales.
	 do uu=1,6
	 do uuu=1,6
     Nn(uu,uuu)=0.0_8
	 enddo
	 do uuu=1,3
	 Nm(uu,uuu)=0.0_8
     Nk(uu,uuu)=0.0_8
     enddo
	 enddo

	 do uu=1,3
	 !Los valores generados en un punto de integraci�n se suman a los generados previamente. As� ser�n calculadas las integrales 
	 !y los valores ir�n en las matrices elementales al salir de este bucle.
	 
	 Mix(1)=((4.0_8*luno(uu)-1.0_8)*yb)/jac
	 Mix(2)=((1.0_8-4.0_8*ldos(uu))*ya)/jac
	 Mix(3)=((4.0_8*luno(uu)+4.0_8*ldos(uu)-3.0_8)*(yb-ya))/jac
	 Mix(4)=(4.0_8*ldos(uu)*yb-4.0_8*luno(uu)*ya)/jac
	 Mix(5)=((8.0_8*ldos(uu)+4.0_8*luno(uu)-4.0_8)*ya-4.0_8*ldos(uu)*yb)/jac
	 Mix(6)=((4.0_8-4.0_8*ldos(uu)-8.0_8*luno(uu))*yb+4.0_8*luno(uu)*ya)/jac
	 
	 Miy(1)=((1.0_8-4.0_8*luno(uu))*xb)/jac
	 Miy(2)=((4.0_8*ldos(uu)-1.0_8)*xa)/jac
	 Miy(3)=((4.0_8*luno(uu)+4.0_8*ldos(uu)-3.0_8)*(xa-xb))/jac
	 Miy(4)=(4.0_8*luno(uu)*xa-4.0_8*ldos(uu)*xb)/jac
	 Miy(5)=(4.0_8*ldos(uu)*xb+(4.0_8-4.0_8*luno(uu)-8.0_8*ldos(uu))*xa)/jac
	 Miy(6)=((8.0_8*luno(uu)+4.0_8*ldos(uu)-4.0_8)*xb-4.0_8*luno(uu)*xa)/jac
	 
	 Mp(1)=luno(uu)
	 Mp(2)=ldos(uu)
	 Mp(3)=1.0_8-luno(uu)-ldos(uu)
	 		
	  do ui=1,6
  	  do uj=1,6
	  !Posicion 1,1 y 2,2. Caja A.
	  Nn(uj,ui)=Nn(uj,ui)+(Mix(uj)*Mix(ui)+Miy(uj)*Miy(ui))*ja(uu)*jac
	  enddo
	  enddo
	  do ui=1,3
  	  do uj=1,6   	 
	  !Posicion 1,3. Caja Bx. 
	  Nm(uj,ui)=Nm(uj,ui)+(Mix(uj)*Mp(ui))*ja(uu)*jac	 				   
	  !Posicion 2,3. Caja By.
	  Nk(uj,ui)=Nk(uj,ui)+(Miy(uj)*Mp(ui))*ja(uu)*jac
	  enddo
	  enddo  
    enddo
	!Escritura de los coeficientes de las matrices elementales en los vectores 'ita,sa' llamando a la subrutina suma.
	!La caja A se escribe en dos posiciones diferentes.
	call suma(nu*Nn/del,n,0,0,3*i,6,6,poi)
	call suma(nu*Nn/del,n,i,i,3*i,6,6,podosi)
	call suma(-Nm*9.81_8/del,n,0,2*i,3*i,6,3,potresi) 
	call suma(-Nk*9.81_8/del,n,i,2*i,3*i,6,3,potresi)
enddo
close(1)
deallocate(poi,podosi,potresi)
end

!----------------------------------------------------------------------------------------------------------------------------------
!Subrutina CAJASBT
!A continuaci�n se calculan coeficientes de las matrices elementales que forman las cajas Bxt y Byt.
!'Nn' son las matrices elementales de la caja Bxt, 'Nm' son las matrices elementales de la caja Byt. Son las transpuestas de las
!cajas Bx y By definidas en la subrutina cajasab.
!Estas cajas aparecen al discretizar la ecuaci�n de continuidad para las ecuaciones Navier-Stokes 2D. Estas ecuaciones pueden ser 
!utilizadas en la primera o primeras iteraciones de Picard para buscar una buena aproximaci�n inicial para las ecuaciones de aguas 
!someras. 
!----------------------------------------------------------------------------------------------------------------------------------
subroutine cajasbt (i,j,x,y,del)
use elemental
integer*4, dimension(:),allocatable::poi,podosi
integer*4 i,j
real*8 x(i),y(i),Nn(3,6),Nm(3,6),Mix(6),Miy(6),Mp(3),del

allocate(poi(3*i),podosi(3*i))

 39     format(4/,A80)
 40     format(6X,6(X,I5))

!Inizializaci�n previa de variables:
!-----------------------------------
luno=(/0.5_8, 0.0_8, 0.5_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8/)		   
ldos=(/0.5_8, 0.5_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8/)
ja=(/1.0_8/6.0_8, 1.0_8/6.0_8, 1.0_8/6.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8/)
do u=1,3*i
poi(u)=posi(u)
podosi(u)=posdosi(u)
enddo

!C�lculo de los valores para cada elemento
!-----------------------------------------
open(unit=1,file='C:\malla.txt',status='old')
read(1,39)ac
 do u=1,j 
     read(1,40) n(1),n(4),n(2),n(5),n(3),n(6)
 	 
	 xa=x(n(1))-x(n(3)) 
	 xb=x(n(2))-x(n(3)) 
	 ya=y(n(1))-y(n(3))
	 yb=y(n(2))-y(n(3))
	 jac=abs(xa*yb-xb*ya)
	   
	 do uu=1,3
     do uuu=1,6
     Nn(uu,uuu)=0.0_8
     Nm(uu,uuu)=0.0_8
     enddo
     enddo

	 do uu=1,3     
	 
	 Mix(1)=((4.0_8*luno(uu)-1.0_8)*yb)/jac
	 Mix(2)=((1.0_8-4.0_8*ldos(uu))*ya)/jac
	 Mix(3)=((4.0_8*luno(uu)+4.0_8*ldos(uu)-3.0_8)*(yb-ya))/jac
	 Mix(4)=(4.0_8*ldos(uu)*yb-4.0_8*luno(uu)*ya)/jac
	 Mix(5)=((8.0_8*ldos(uu)+4.0_8*luno(uu)-4.0_8)*ya-4.0_8*ldos(uu)*yb)/jac
	 Mix(6)=((4.0_8-4.0_8*ldos(uu)-8.0_8*luno(uu))*yb+4.0_8*luno(uu)*ya)/jac
	
	 Miy(1)=((1.0_8-4.0_8*luno(uu))*xb)/jac
	 Miy(2)=((4.0_8*ldos(uu)-1.0_8)*xa)/jac
	 Miy(3)=((4.0_8*luno(uu)+4.0_8*ldos(uu)-3.0_8)*(xa-xb))/jac
	 Miy(4)=(4.0_8*luno(uu)*xa-4.0_8*ldos(uu)*xb)/jac
	 Miy(5)=(4.0_8*ldos(uu)*xb+(4.0_8-4.0_8*luno(uu)-8.0_8*ldos(uu))*xa)/jac
	 Miy(6)=((8.0_8*luno(uu)+4.0_8*ldos(uu)-4.0_8)*xb-4.0_8*luno(uu)*xa)/jac
	 
	 Mp(1)=luno(uu)
	 Mp(2)=ldos(uu)
	 Mp(3)=1.0_8-luno(uu)-ldos(uu)
	 
	  do ui=1,6
  	  do uj=1,3   	 
	  !Posicion 3,1. Caja Bxt.	  
	  Nn(uj,ui)=Nn(uj,ui)+(Mix(ui)*Mp(uj))*ja(uu)*jac 	   
	  !Posicion 3,2. Caja Byt.
	  Nm(uj,ui)=Nm(uj,ui)+(Miy(ui)*Mp(uj))*ja(uu)*jac
	  enddo
	  enddo	   	  
     enddo 
	 call suma(Nn/del,n,2*i,0,3*i,3,6,poi)
	 call suma(Nm/del,n,2*i,i,3*i,3,6,podosi)
 enddo
 close(1)
 deallocate(poi,podosi) 
end

!--------------------------------------------------------------------------------------------------------------------------------------
!Subrutina CAJASDE
!A continuaci�n se calculan coeficientes de las matrices elementales que forman las cajas Dx y Dy.
!'Nn' son las matrices elementales de la caja Dx, 'Nm' son las matrices elementales de la caja Dy. Estas cajas aparecen al discretizar 
!la ecuaci�n de continuidad para las ecuaciones de aguas someras. 
!En 'vt' est� el calado soluci�n de la iteraci�n anterior del que dependen estas cajas. Por ello se trata de una caja no lineal del 
!sistema. 'hxy' es el calado evaluado en cualquier punto de integraci�n del elemento, a partir de los valores de calado de la 
!iteraci�n anterior. Aqu� s�lo se pueden utilizar valores positivos de calado. 
!--------------------------------------------------------------------------------------------------------------------------------------
subroutine cajasde(i,j,x,y,vt,del)
use elemental
integer*4, dimension(:),allocatable::poi,podosi
integer*4 i,j										   
real*8 x(i),y(i),Nn(3,6),Nm(3,6),Mi(6),Mp(3),Mpx(3),Mpy(3),hxy,vt(3*i),del

allocate(poi(3*i),podosi(3*i))
 
39     format(4/,A80)
40     format(6X,6(X,I5))

!Inizializaci�n previa de variables:
!-----------------------------------
luno=(/1.0_8/3.0_8, 0.6_8, 0.2_8, 0.2_8, 0.0_8, 0.0_8, 0.0_8/)		   
ldos=(/1.0_8/3.0_8, 0.2_8, 0.6_8, 0.2_8, 0.0_8, 0.0_8, 0.0_8/)
ja=(/-27.0_8/96.0_8, 25.0_8/96.0_8, 25.0_8/96.0_8, 25.0_8/96.0_8, 0.0_8, 0.0_8, 0.0_8/)
do u=1,3*i
poi(u)=posi(u)
podosi(u)=posdosi(u)
enddo

!C�lculo de los valores para cada elemento
!-----------------------------------------
open(unit=1,file='C:\malla.txt',status='old')
read(1,39)ac
 do u=1,j
     read(1,40) n(1),n(4),n(2),n(5),n(3),n(6)
  
     xa=x(n(1))-x(n(3)) 
	 xb=x(n(2))-x(n(3)) 
	 ya=y(n(1))-y(n(3))
	 yb=y(n(2))-y(n(3))
	 jac=abs(xa*yb-xb*ya)

	 do uu=1,3
     do uuu=1,6
     Nn(uu,uuu)=0.0_8
	 Nm(uu,uuu)=0.0_8
     enddo
     enddo

	 !Integraci�n con cuatro puntos. Polinomios de grado superior a grado 3 quedar�an integrados con error.
	 do uu=1,4
	 			   
	 Mi(1)=luno(uu)*(2.0_8*luno(uu)-1.0_8)
	 Mi(2)=ldos(uu)*(2.0_8*ldos(uu)-1.0_8)
	 Mi(3)=(1.0_8-luno(uu)-ldos(uu))*(1.0_8-2.0_8*luno(uu)-2.0_8*ldos(uu))
	 Mi(4)=4.0_8*luno(uu)*ldos(uu)
	 Mi(5)=4.0_8*ldos(uu)*(1.0_8-luno(uu)-ldos(uu))
	 Mi(6)=4.0_8*luno(uu)*(1.0_8-luno(uu)-ldos(uu))

	 Mp(1)=luno(uu)
	 Mp(2)=ldos(uu)
	 Mp(3)=1.0_8-luno(uu)-ldos(uu)

	 Mpx(1)=(y(n(2))-y(n(3)))/jac
	 Mpx(2)=(y(n(3))-y(n(1)))/jac
	 Mpx(3)=(y(n(1))-y(n(2)))/jac
	 
	 Mpy(1)=(x(n(3))-x(n(2)))/jac
	 Mpy(2)=(x(n(1))-x(n(3)))/jac
	 Mpy(3)=(x(n(2))-x(n(1)))/jac 
	 
	 !Se utilizan los calados calculados en la iteracion anterior.
	 hxy=vt(2*i+n(1))*Mp(1)+vt(2*i+n(2))*Mp(2)+vt(2*i+n(3))*Mp(3)
	
	 do ui=1,6
	 do uj=1,3	 
	 !Posicion 3,1. Caja Dx.
	 Nn(uj,ui)=Nn(uj,ui)+(Mpx(uj)*Mi(ui)*hxy)*ja(uu)*jac	 
	 !Posicion 3,2. Caja Dy.
	 Nm(uj,ui)=Nm(uj,ui)+(Mpy(uj)*Mi(ui)*hxy)*ja(uu)*jac
	 enddo
	 enddo		   	  
   enddo
   call suma(-Nn/del,n,2*i,0,3*i,3,6,poi)
   call suma(-Nm/del,n,2*i,i,3*i,3,6,podosi)
 enddo
 close(1)
 deallocate(poi,podosi)
end

!------------------------------------------------------------------------------------------------------------------------------------------
!Subrutina MATRIZNOLINEAL
!A continuaci�n se calculan coeficientes de las matrices elementales que forman la caja C.
!'Nn' son las matrices elementales de esta caja. Estas cajas aparecen al discretizar las ecuaciones din�micas (ec. para flujo superficial).
!En 'vv' est�n las velocidades soluci�n de la iteraci�n anterior de las que dependen estas cajas. Por ello se trata de una caja no 
!lineal del sistema.
!'Mu, Mv' son las velocidades en direcci�n x y en direcci�n y evaluadas en cualquier punto de integraci�n del elemento, a partir de 
!los valores de velocidad la iteraci�n anterior.
!------------------------------------------------------------------------------------------------------------------------------------------
subroutine matriznolineal(i,j,x,y,vv,del)
use elemental
integer*4, dimension(:),allocatable::poi,podosi
integer*4 i,j
real*8 x(i),y(i),Nn(6,6),Mi(6),Mix(6),Miy(6),Mu,Mv,vv(3*i),del

allocate(poi(3*i),podosi(3*i))

 39     format(4/,A80)
 40     format(6X,6(X,I5))

!Inizializaci�n previa de variables:
!-----------------------------------
luno=(/1.0_8/3.0_8, 0.0597158717_8, 0.4701420641_8, 0.4701420641_8, 0.7974269853_8, 0.1012865073_8, 0.1012865073_8/)		   
ldos=(/1.0_8/3.0_8, 0.4701420641_8, 0.0597158717_8, 0.4701420641_8, 0.1012865073_8, 0.7974269853_8, 0.1012865073_8/)
ja=(/0.1125_8, 0.06619707635_8, 0.06619707635_8, 0.06619707635_8, 0.06296959025_8, 0.06296959025_8, 0.06296959025_8/)  
do u=1,3*i
poi(u)=posi(u)
podosi(u)=posdosi(u)
enddo

!C�lculo de los valores para cada elemento
!-----------------------------------------
open(unit=1,file='C:\malla.txt',status='old')
read(1,39)ac
 do u=1,j
     read(1,40) n(1),n(4),n(2),n(5),n(3),n(6)

	 xa=x(n(1))-x(n(3)) 
	 xb=x(n(2))-x(n(3)) 
	 ya=y(n(1))-y(n(3))
	 yb=y(n(2))-y(n(3))  
	 jac=abs(xa*yb-xb*ya)
	 	 
	 do uu=1,6
     do uuu=1,6
     Nn(uu,uuu)=0.0_8
     enddo
     enddo

	 !Integraci�n gauss de cuarto orden con siete puntos de integraci�n. Polinomios de grado superior a grado 5 quedar�an integrados con error.
	 do uu=1,7
	 Mi(1)=luno(uu)*(2.0_8*luno(uu)-1.0_8)
	 Mi(2)=ldos(uu)*(2.0_8*ldos(uu)-1.0_8)
	 Mi(3)=(1.0_8-luno(uu)-ldos(uu))*(1.0_8-2.0_8*luno(uu)-2.0_8*ldos(uu))
	 Mi(4)=4.0_8*luno(uu)*ldos(uu)
	 Mi(5)=4.0_8*ldos(uu)*(1.0_8-luno(uu)-ldos(uu))
	 Mi(6)=4.0_8*luno(uu)*(1.0_8-luno(uu)-ldos(uu))

	 Mix(1)=((4.0_8*luno(uu)-1.0_8)*yb)/jac
	 Mix(2)=((1.0_8-4.0_8*ldos(uu))*ya)/jac
	 Mix(3)=((4.0_8*luno(uu)+4.0_8*ldos(uu)-3.0_8)*(yb-ya))/jac
	 Mix(4)=(4.0_8*ldos(uu)*yb-4.0_8*luno(uu)*ya)/jac
	 Mix(5)=((8.0_8*ldos(uu)+4.0_8*luno(uu)-4.0_8)*ya-4.0_8*ldos(uu)*yb)/jac
	 Mix(6)=((4.0_8-4.0_8*ldos(uu)-8.0_8*luno(uu))*yb+4.0_8*luno(uu)*ya)/jac
	
	 Miy(1)=((1.0_8-4.0_8*luno(uu))*xb)/jac
	 Miy(2)=((4.0_8*ldos(uu)-1.0_8)*xa)/jac
	 Miy(3)=((4.0_8*luno(uu)+4.0_8*ldos(uu)-3.0_8)*(xa-xb))/jac
	 Miy(4)=(4.0_8*luno(uu)*xa-4.0_8*ldos(uu)*xb)/jac
	 Miy(5)=(4.0_8*ldos(uu)*xb+(4.0_8-4.0_8*luno(uu)-8.0_8*ldos(uu))*xa)/jac
	 Miy(6)=((8.0_8*luno(uu)+4.0_8*ldos(uu)-4.0_8)*xb-4.0_8*luno(uu)*xa)/jac

	 !Se calcula el producto escalar
	 Mu=dot_product(vv(n),Mi)   
	 Mv=dot_product(vv(n+i),Mi) 	 	
	 
	  do ui=1,6
  	  do uj=1,6 
	  !Posici�n 1,1 y 2,2. Caja C.
	  Nn(uj,ui)=Nn(uj,ui)+Mi(uj)*(Mix(ui)*Mu+Miy(ui)*Mv)*ja(uu)*jac
	  enddo
	  enddo	 
	 enddo
	 !La caja C se escribe en dos posiciones diferentes. Adem�s se trata de las posiciones donde se escribi� la caja A en la 
	 !subrutina cajasab. De ah� la utilidad de las variables 'poi,podosi'.
	 call suma(Nn/del,n,0,0,3*i,6,6,poi)
	 call suma(Nn/del,n,i,i,3*i,6,6,podosi)
 enddo
 close(1)
 deallocate(poi,podosi)
end
	
!------------------------------------------------------------------------------------------------------------------------------------------
!Subrutina TIMEASU
!A continuaci�n se calculan coeficientes de las matrices elementales que forman la caja M, t�picamente llamada matriz de masa.
!'Nn' son las matrices elementales de esta caja. Estas cajas aparecen al discretizar las ecuaciones din�micas (ec. para flujo superficial).
!S�lo entra aqu� en caso de resolver el problema no estacionario. 
!------------------------------------------------------------------------------------------------------------------------------------------
subroutine timeasu(i,j,x,y,At)
use elemental
integer*4, dimension(:),allocatable::poi,podosi
integer*4 i,j				  
real*8 x(i),y(i),Nn(6,6),Mi(6),At

allocate(poi(3*i),podosi(3*i))

 39     format(4/,A80)
 40     format(6X,6(X,I5))

!Inizializaci�n previa de variables:
!-----------------------------------
luno=(/1.0_8/3.0_8, 0.0597158717_8, 0.4701420641_8, 0.4701420641_8, 0.7974269853_8, 0.1012865073_8, 0.1012865073_8/)		   
ldos=(/1.0_8/3.0_8, 0.4701420641_8, 0.0597158717_8, 0.4701420641_8, 0.1012865073_8, 0.7974269853_8, 0.1012865073_8/)
ja=(/0.1125_8, 0.06619707635_8, 0.06619707635_8, 0.06619707635_8, 0.06296959025_8, 0.06296959025_8, 0.06296959025_8/)
do u=1,3*i
poi(u)=posi(u)
podosi(u)=posdosi(u)
enddo

!C�lculo de los valores para cada elemento
!-----------------------------------------
open(unit=1,file='C:\malla.txt',status='old')  
read(1,39)ac	 
 do u=1,j
     read(1,40) n(1),n(4),n(2),n(5),n(3),n(6)
  
	 xa=x(n(1))-x(n(3)) 
	 xb=x(n(2))-x(n(3)) 
	 ya=y(n(1))-y(n(3))
	 yb=y(n(2))-y(n(3))  	 
	 jac=abs(xa*yb-xb*ya)
	 	 
	 do uu=1,6
     do uuu=1,6
     Nn(uu,uuu)=0.0_8
     enddo
     enddo

	 do uu=1,7
	 Mi(1)=luno(uu)*(2.0_8*luno(uu)-1.0_8)
	 Mi(2)=ldos(uu)*(2.0_8*ldos(uu)-1.0_8)
	 Mi(3)=(1.0_8-luno(uu)-ldos(uu))*(1.0_8-2.0_8*luno(uu)-2.0_8*ldos(uu))
	 Mi(4)=4.0_8*luno(uu)*ldos(uu)
	 Mi(5)=4.0_8*ldos(uu)*(1.0_8-luno(uu)-ldos(uu))
	 Mi(6)=4.0_8*luno(uu)*(1.0_8-luno(uu)-ldos(uu))

	  do ui=1,6
  	  do uj=1,6 
	  !Posici�n 1,1 y 2,2. Caja M.
	  Nn(uj,ui)=Nn(uj,ui)+(Mi(uj)*Mi(ui))*ja(uu)*jac
	  enddo
	  enddo	 
	 enddo
	 call suma(Nn/At,n,0,0,3*i,6,6,poi)
	 call suma(Nn/At,n,i,i,3*i,6,6,podosi)
 enddo
 close(1)
 deallocate(poi,podosi)
end

!------------------------------------------------------------------------------------------------------------------------------
!Subrutina TIMEASD
!A continuaci�n se calculan coeficientes de las matrices elementales que forman la caja N, t�picamente llamada matriz de masa.
!'Nm' son las matrices elementales de esta caja. Estas cajas aparecen al discretizar la ecuaci�n de continuidad de las 
!ecuaciones de aguas someras.
!S�lo entra aqu� en caso de resolver el problema no estacionario.
!------------------------------------------------------------------------------------------------------------------------------
subroutine timeasd(i,j,x,y,At)
use elemental
integer*4, dimension(:),allocatable::potresi
integer*4 i,j				  
real*8 x(i),y(i),Nm(3,3),Mp(3),At

allocate(potresi(3*i))

 39     format(4/,A80)
 40     format(6X,6(X,I5))

!Inizializaci�n previa de variables:
!-----------------------------------
luno=(/0.5_8, 0.0_8, 0.5_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8/)		   
ldos=(/0.5_8, 0.5_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8/)
ja=(/1.0_8/6.0_8, 1.0_8/6.0_8, 1.0_8/6.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8/)
do u=1,3*i
potresi(u)=postresi(u)
enddo

!C�lculo de los valores para cada elemento
!-----------------------------------------
open(unit=1,file='C:\malla.txt',status='old')  
read(1,39)ac	 
 do u=1,j
     read(1,40) n(1),n(4),n(2),n(5),n(3),n(6)
  
	 xa=x(n(1))-x(n(3)) 
	 xb=x(n(2))-x(n(3)) 
	 ya=y(n(1))-y(n(3))
	 yb=y(n(2))-y(n(3))  	 
	 jac=abs(xa*yb-xb*ya)
	 	 
	 do uu=1,3
     do uuu=1,3
     Nm(uu,uuu)=0.0_8
     enddo
     enddo

	 do uu=1,3
	 Mp(1)=luno(uu)
	 Mp(2)=ldos(uu)
	 Mp(3)=1.0_8-luno(uu)-ldos(uu)

	  do ui=1,3
  	  do uj=1,3
      !Posici�n 3,3. Caja N.
	  Nm(uj,ui)=Nm(uj,ui)+(Mp(uj)*Mp(ui))*ja(uu)*jac
	  enddo
	  enddo	 
	 enddo
	 call suma(Nm/At,n,2*i,2*i,3*i,3,3,potresi)
 enddo
 close(1)
 deallocate(potresi)
end

!-----------------------------------------------------------------------------------------------------------------------------------------
!Subrutina F
!Con la siguiente subrutina se calcula la fricci�n por la tensi�n de fondo a trav�s del n�mero de Manning. 
!Se calculan vectores elementales que ir�n al t�rmino independiente del sistema y que aparecen en las ecuaciones din�micas de las 
!ecuaciones de aguas someras.
!Con 'yn' se da un rozamiento tal que I=i. As�, independientemente de la velocidad y el calado, se dar� un calado normal.
!Esto es equivalente a resolver un sistema con la ecuaci�n de continuidad de aguas someras y las ecuaciones din�micas de Navier-Stokes 2D.
!-----------------------------------------------------------------------------------------------------------------------------------------
subroutine f(i,j,x,y,z,ma,yn,vv,Nx,Ny)	 
use elemental
integer*4 i,j	 
real*8 x(i),y(i),Nx(i),Ny(i),s,maning,ma(i),vv(3*i),Mi(6),Mu,Mv,Mui,Mvi,Mhi,dist,a,b,c
real*8 z(i),Mix(6),Miy(6),Mpx,Mpy  
character yn*2
			 
39     format(4/,A80)
40     format(6X,6(X,I5))
		  													 
!Inizializaci�n previa de variables:
!-----------------------------------
luno=(/1.0_8/3.0_8, 0.0597158717_8, 0.4701420641_8, 0.4701420641_8, 0.7974269853_8, 0.1012865073_8, 0.1012865073_8/)		   
ldos=(/1.0_8/3.0_8, 0.4701420641_8, 0.0597158717_8, 0.4701420641_8, 0.1012865073_8, 0.7974269853_8, 0.1012865073_8/)
ja=(/0.1125_8, 0.06619707635_8, 0.06619707635_8, 0.06619707635_8, 0.06296959025_8, 0.06296959025_8, 0.06296959025_8/)    		   

!C�lculo de los valores para cada elemento
!-----------------------------------------
open(unit=1,file='C:\malla.txt',status='old')
read(1,39)ac
 do u=1,j
 read(1,40) n(1),n(4),n(2),n(5),n(3),n(6)

	 xa=x(n(1))-x(n(3)) 
	 xb=x(n(2))-x(n(3)) 
	 ya=y(n(1))-y(n(3))
	 yb=y(n(2))-y(n(3))  
	 jac=abs(xa*yb-xb*ya)	

	 !C�lculo del m�dulo de la velocidad en cada nodo esquina.
	 a=sqrt(vv(n(1))**2+vv(i+n(1))**2)
	 b=sqrt(vv(n(2))**2+vv(i+n(2))**2)
	 c=sqrt(vv(n(3))**2+vv(i+n(3))**2)
	 !La distancia elemental se estima como el di�metro del c�rculo de �rea igual a la del elemento.
	 dist=sqrt((2.0_8*jac)/(3.14159_8))
	 
	 !Valores medios de velocidad, calado y manning en el elemento (valores en el centro del elemento).
	 Mui=-(vv(n(1))+vv(n(2))+vv(n(3)))/9.0_8+(vv(n(4))+vv(n(5))+vv(n(6)))*4.0_8/9.0_8   
	 Mvi=-(vv(n(1)+i)+vv(n(2)+i)+vv(n(3)+i))/9.0_8+(vv(n(4)+i)+vv(n(5)+i)+vv(n(6)+i))*4.0_8/9.0_8
	 Mhi=(vv(n(1)+2*i)+vv(n(2)+2*i)+vv(n(3)+2*i))/3.0_8
	 maning=(ma(n(1))+ma(n(2))+ma(n(3)))/3.0_8
	 
	 !Se da un coeficiente de Manning atendiendo al tipo de elemento.
	 if (((a.eq.0.0).and.(b.eq.0.0)).or.((a.eq.0.0).and.(c.eq.0.0)).or.((b.eq.0.0).and.(c.eq.0.0))) then
	  !Elementos pegados a un contorno con velocidades nulas. Se consideran s�lo elementos con dos nodos apoyados en la pared. 
	  !El per�metro mojado tiene en cuenta la pared lateral y el radio hidr�ulico es funci�n de la distancia. 
	  s=((maning**2.0_8))/(((dist*Mhi)/(dist+Mhi))**(4.0_8/3.0_8))	  
	 else 
	  !Se consideran elementos con un nodo apoyado en la pared o con ning�n nodo apoyado.
	  !En aqu�llos con un nodo apoyado se tiene la pared en un punto infinitesimal del elemento. �sta ya se considera en elementos con dos 
	  !nodos apoyados. El radio hidr�ulico se estima como el calado medio.
	  s=((maning**2.0_8))/(Mhi**(4.0_8/3.0_8)) 
	 endif
	  
	 do uu=1,7	      
	 Mi(1)=luno(uu)*(2.0_8*luno(uu)-1.0_8)
	 Mi(2)=ldos(uu)*(2.0_8*ldos(uu)-1.0_8)
	 Mi(3)=(1.0_8-luno(uu)-ldos(uu))*(1.0_8-2.0_8*luno(uu)-2.0_8*ldos(uu))
	 Mi(4)=4.0_8*luno(uu)*ldos(uu)
	 Mi(5)=4.0_8*ldos(uu)*(1.0_8-luno(uu)-ldos(uu))
	 Mi(6)=4.0_8*luno(uu)*(1.0_8-luno(uu)-ldos(uu))

	 !Procedimiento con poco sentido f�sico.
	 !Se obliga a que Ix=ix, Iy=iy (con i cambiado de signo, ix ser� positivo si z baja a medida que crecen los valores de x).
	 !Se ha comprobado que efectivamente con esto la altura de la l�mina de agua se ajusta a la cota del terreno (calado normal).
	 if (yn.eq.'si') then
	 Mix(1)=((4.0_8*luno(uu)-1.0_8)*yb)/jac
	 Mix(2)=((1.0_8-4.0_8*ldos(uu))*ya)/jac
	 Mix(3)=((4.0_8*luno(uu)+4.0_8*ldos(uu)-3.0_8)*(yb-ya))/jac
	 Mix(4)=(4.0_8*ldos(uu)*yb-4.0_8*luno(uu)*ya)/jac
	 Mix(5)=((8.0_8*ldos(uu)+4.0_8*luno(uu)-4.0_8)*ya-4.0_8*ldos(uu)*yb)/jac
	 Mix(6)=((4.0_8-4.0_8*ldos(uu)-8.0_8*luno(uu))*yb+4.0_8*luno(uu)*ya)/jac
	
	 Miy(1)=((1.0_8-4.0_8*luno(uu))*xb)/jac
	 Miy(2)=((4.0_8*ldos(uu)-1.0_8)*xa)/jac
	 Miy(3)=((4.0_8*luno(uu)+4.0_8*ldos(uu)-3.0_8)*(xa-xb))/jac
	 Miy(4)=(4.0_8*luno(uu)*xa-4.0_8*ldos(uu)*xb)/jac
	 Miy(5)=(4.0_8*ldos(uu)*xb+(4.0_8-4.0_8*luno(uu)-8.0_8*ldos(uu))*xa)/jac
	 Miy(6)=((8.0_8*luno(uu)+4.0_8*ldos(uu)-4.0_8)*xb-4.0_8*luno(uu)*xa)/jac

	  Mpx=dot_product(z(n),Mix)	  
	  Mpy=dot_product(z(n),Miy)
	  do uuu=1,6
	  Nx(n(uuu))=Nx(n(uuu))-Mi(uuu)*Mpx*jac*ja(uu) 	 
	  Ny(n(uuu))=Ny(n(uuu))-Mi(uuu)*Mpy*jac*ja(uu) 
	  enddo
	 
	 !Procedimiento para calcular con las ecuaciones de aguas someras. 
	 !Se calcula la pendiente de fricci�n a partir de las variables de velocidad y calado.
	 else
	  Mu=dot_product(vv(n),Mi)   
	  Mv=dot_product(vv(n+i),Mi) 	  
	  do uuu=1,6
	  Nx(n(uuu))=Nx(n(uuu))+(Mi(uuu)*Mu*sqrt(Mui**2+Mvi**2)*s)*ja(uu)*jac 	 
	  Ny(n(uuu))=Ny(n(uuu))+(Mi(uuu)*Mv*sqrt(Mui**2+Mvi**2)*s)*ja(uu)*jac 
	  enddo	
	 endif
	 enddo
 enddo 
 close(1)
 end

!------------------------------------------------------------------------------------------------------------------------
!Subrutina JACOB
!Esta subrutina calcula los t�rminos que faltan de la matriz jacobiana del sistema para las ecuaciones de aguas someras.
!Tambi�n se dan condiciones de contorno sobre el sistema como en la subrutina reducciondelsistema.
!As�, se tendr� preparado el sistema a resolver para el m�todo de Newton.
!------------------------------------------------------------------------------------------------------------------------
subroutine jacob (i,j,x,y,vv,nu,man,vector,v,vic,vb,cc,ten,bcg,ndim)
use elemental
use allocatacion
integer*4, dimension(:),allocatable::poi,podosi,potresi,ck
integer*4 i,j,ndim,v(i),cc,inc,w,k
real*8 x(i),y(i),Ncua(6,6),Ncub(6,6),Ncva(6,6),Ncvb(6,6),Nbhx(6,3),Nbhy(6,3),Nbuh(3,6),Nbvh(3,6),Nbu(3,3),Nbv(3,3)
real*8 Ndyx(3,3),Nfx(6,6),Nfy(6,6),Nbnux(6,6),Mi(6),Mix(6),Miy(6),Mp(3),Mpx(3),Mpy(3),vv(3*i),vector(3*i),vic(3*i)
real*8 Mu,Mv,Mux,Mvx,Muy,Mvy,Mui,Mvi,Mhi,maning,man(i),amn,bmn,akm,bkm,ank,bnk,a,b,c,dist,s,ma,mb,mc,md,me,nu,vb(3*i)
character ten*2,bcg*2

allocate(poi(3*i),podosi(3*i),potresi(3*i),ck(3*i))

 39     format(4/,A80)										  
 40     format(6X,6(X,I5))

!Inizializaci�n previa de variables:
!-----------------------------------
!Se utilizar�n dos integraciones diferentes en el dominio, una con 'luno,ldos,ja' y otra con 'lunot,ldost,jat'.
luno=(/1.0_8/3.0_8, 0.6_8, 0.2_8, 0.2_8, 0.0_8, 0.0_8, 0.0_8/)		   
ldos=(/1.0_8/3.0_8, 0.2_8, 0.6_8, 0.2_8, 0.0_8, 0.0_8, 0.0_8/)
ja=(/-27.0_8/96.0_8, 25.0_8/96.0_8, 25.0_8/96.0_8, 25.0_8/96.0_8, 0.0_8, 0.0_8, 0.0_8/)
lunot=(/1.0_8/3.0_8, 0.0597158717_8, 0.4701420641_8, 0.4701420641_8, 0.7974269853_8, 0.1012865073_8, 0.1012865073_8/)		   
ldost=(/1.0_8/3.0_8, 0.4701420641_8, 0.0597158717_8, 0.4701420641_8, 0.1012865073_8, 0.7974269853_8, 0.1012865073_8/)
jat=(/0.1125_8, 0.06619707635_8, 0.06619707635_8, 0.06619707635_8, 0.06296959025_8, 0.06296959025_8, 0.06296959025_8/)
do u=1,3*i
poi(u)=posi(u)
podosi(u)=posdosi(u)
potresi(u)=postresi(u)
enddo
ma=9.0_8
mb=12.0_8
mc=8.0_8
md=3.0_8
me=4.0_8

open(unit=1,file='C:\malla.txt',status='old')
read(1,39)ac
 do u=1,j
     read(1,40) n(1),n(4),n(2),n(5),n(3),n(6)

	 xa=x(n(1))-x(n(3)) 
	 xb=x(n(2))-x(n(3)) 
	 ya=y(n(1))-y(n(3))
	 yb=y(n(2))-y(n(3))  
	 jac=abs(xa*yb-xb*ya)
	 
	 !C�lculo de par�metros para el c�lculo del t�rmino correspondiente a la pendiente de fricci�n, de un modo similar al de la subrutina f.
	 !Se considerar� s�lo el caso yn='no' indicado en la subrutina f.
	 a=sqrt(vv(n(1))**2+vv(i+n(1))**2)
	 b=sqrt(vv(n(2))**2+vv(i+n(2))**2)
	 c=sqrt(vv(n(3))**2+vv(i+n(3))**2)
	 dist=sqrt((2.0_8*jac)/(3.14159_8))
	 Mui=-(vv(n(1))+vv(n(2))+vv(n(3)))/9.0_8+(vv(n(4))+vv(n(5))+vv(n(6)))*4.0_8/9.0_8   
	 Mvi=-(vv(n(1)+i)+vv(n(2)+i)+vv(n(3)+i))/9.0_8+(vv(n(4)+i)+vv(n(5)+i)+vv(n(6)+i))*4.0_8/9.0_8
	 Mhi=(vv(n(1)+2*i)+vv(n(2)+2*i)+vv(n(3)+2*i))/3.0_8
	 maning=(man(n(1))+man(n(2))+man(n(3)))/3.0_8
	  if (((a.eq.0.0).and.(b.eq.0.0)).or.((a.eq.0.0).and.(c.eq.0.0)).or.((b.eq.0.0).and.(c.eq.0.0))) then
	   s=((maning**2.0_8))/(((dist*Mhi)/(dist+Mhi))**(4.0_8/3.0_8))	  
	  else 
	   s=((maning**2.0_8))/(Mhi**(4.0_8/3.0_8)) 
	  endif

	 do uu=1,6
     do uuu=1,6
	 Ncua(uu,uuu)=0.0_8
	 Ncub(uu,uuu)=0.0_8
	 Ncva(uu,uuu)=0.0_8
	 Ncvb(uu,uuu)=0.0_8 
	 Nfx(uu,uuu)=0.0_8
	 Nfy(uu,uuu)=0.0_8	
     enddo
     enddo
	 do uu=1,3
     do uuu=1,3
	 Ndyx(uu,uuu)=0.0_8	
     enddo
     enddo

	 call direccion(x,y,i,n(2),n(3),n(1),amn,bmn)
	 call direccion(x,y,i,n(1),n(2),n(3),akm,bkm)
	 call direccion(x,y,i,n(3),n(1),n(2),ank,bnk)
	 
	 !Cajas formadas a partir de integrales de contorno
	 !-------------------------------------------------
	 !Derivadas de las integrales de contorno de los t�rminos de presi�n para el m�todo de Newton.
	 Nbhx(1,1)=(akm*abs(ya-yb)+ank*abs(ya))/6.0_8	 
	 Nbhx(2,2)=(akm*abs(ya-yb)+amn*abs(yb))/6.0_8	 
	 Nbhx(3,3)=(amn*abs(yb)+ank*abs(ya))/6.0_8
	 Nbhx(4,1)=akm*abs(ya-yb)/3.0_8
	 Nbhx(4,2)=akm*abs(ya-yb)/3.0_8
	 Nbhx(5,2)=amn*abs(yb)/3.0_8
	 Nbhx(5,3)=amn*abs(yb)/3.0_8
	 Nbhx(6,1)=ank*abs(ya)/3.0_8
	 Nbhx(6,3)=ank*abs(ya)/3.0_8
	 
	 Nbhy(1,1)=(bkm*abs(xa-xb)+bnk*abs(xa))/6.0_8	 
	 Nbhy(2,2)=(bkm*abs(xa-xb)+bmn*abs(xb))/6.0_8	 
	 Nbhy(3,3)=(bmn*abs(xb)+bnk*abs(xa))/6.0_8
	 Nbhy(4,1)=bkm*abs(xa-xb)/3.0_8
	 Nbhy(4,2)=bkm*abs(xa-xb)/3.0_8
	 Nbhy(5,2)=bmn*abs(xb)/3.0_8
	 Nbhy(5,3)=bmn*abs(xb)/3.0_8
	 Nbhy(6,1)=bnk*abs(xa)/3.0_8
	 Nbhy(6,3)=bnk*abs(xa)/3.0_8

	 !Derivadas de las integrales de contorno de los t�rminos viscosos para el m�todo de Newton. 
	 if (ten.eq.'si') then 
	 !Se consideran si tambi�n se consideran estas integrales de contorno
	 Nbnux(1,1)=md*(yb*ank*abs(ya)+ya*akm*abs(ya-yb)-xb*bnk*abs(xa)-xa*bkm*abs(xa-xb))/jac
	 Nbnux(1,2)=(ya*akm*abs(ya-yb)-xa*bkm*abs(xa-xb))/jac
	 Nbnux(1,3)=(yb*ank*abs(ya)-xb*bnk*abs(xa))/jac
     Nbnux(1,4)=me*(-ya*akm*abs(ya-yb)+xa*bkm*abs(xa-xb))/jac
     Nbnux(1,6)=me*(-yb*ank*abs(ya)+xb*bnk*abs(xa))/jac
	 Nbnux(2,1)=(-ya*akm*abs(ya-yb)+xa*bkm*abs(xa-xb))/jac
	 Nbnux(2,2)=md*(-ya*akm*abs(ya-yb)-(ya-yb)*amn*abs(yb)+xa*bkm*abs(xa-xb)+(xa-xb)*bmn*abs(xb))/jac
	 Nbnux(2,3)=(-(ya-yb)*amn*abs(yb)+(xa-xb)*bmn*abs(xb))/jac
	 Nbnux(2,4)=me*(ya*akm*abs(ya-yb)-xa*bkm*abs(xa-xb))/jac
	 Nbnux(2,5)=me*((ya-yb)*amn*abs(yb)-(xa-xb)*bmn*abs(xb))/jac
	 Nbnux(3,1)=(-yb*ank*abs(ya)+xb*bnk*abs(xa))/jac
	 Nbnux(3,2)=((ya-yb)*amn*abs(yb)-(xa-xb)*bmn*abs(xb))/jac
	 Nbnux(3,3)=md*(-yb*ank*abs(ya)+(ya-yb)*amn*abs(yb)+xb*bnk*abs(xa)-(xa-xb)*bmn*abs(xb))/jac
	 Nbnux(3,5)=me*(-(ya-yb)*amn*abs(yb)+(xa-xb)*bmn*abs(xb))/jac
	 Nbnux(3,6)=me*(yb*ank*abs(ya)-xb*bnk*abs(xa))/jac
     Nbnux(4,1)=me*(ya*akm*abs(ya-yb)-xa*bkm*abs(xa-xb))/jac
     Nbnux(4,2)=me*(-ya*akm*abs(ya-yb)+xa*bkm*abs(xa-xb))/jac
	 Nbnux(5,2)=me*(-(ya-yb)*amn*abs(yb)+(xa-xb)*bmn*abs(xb))/jac
	 Nbnux(5,3)=me*((ya-yb)*amn*abs(yb)-(xa-xb)*bmn*abs(xb))/jac
	 Nbnux(6,1)=me*(yb*ank*abs(ya)-xb*bnk*abs(xa))/jac
	 Nbnux(6,3)=me*(-yb*ank*abs(ya)+xb*bnk*abs(xa))/jac
	 else
	 do uu=1,6
     do uuu=1,6
	 Nbnux(uu,uuu)=0.0_8	
     enddo
     enddo
	 endif

	 !Derivadas de las integrales de contorno de la ecuaci�n de continuidad para el m�todo de Newton.
	 Nbuh(1,1)=(ma*vv(2*i+n(1))+vv(2*i+n(2)))*akm*abs(ya-yb)+(ma*vv(2*i+n(1))+vv(2*i+n(3)))*ank*abs(ya)
	 Nbuh(1,2)=(-vv(2*i+n(1))+vv(2*i+n(2)))*akm*abs(ya-yb)
	 Nbuh(1,3)=(-vv(2*i+n(1))+vv(2*i+n(3)))*ank*abs(ya)
	 Nbuh(1,4)=(mb*vv(2*i+n(1))+mc*vv(2*i+n(2)))*akm*abs(ya-yb)
	 Nbuh(1,6)=(mb*vv(2*i+n(1))+mc*vv(2*i+n(3)))*ank*abs(ya)
	 Nbuh(2,1)=(vv(2*i+n(1))-vv(2*i+n(2)))*akm*abs(ya-yb)
	 Nbuh(2,2)=(vv(2*i+n(1))+ma*vv(2*i+n(2)))*akm*abs(ya-yb)+(ma*vv(2*i+n(2))+vv(2*i+n(3)))*amn*abs(yb)
	 Nbuh(2,3)=(-vv(2*i+n(2))+vv(2*i+n(3)))*amn*abs(yb)
	 Nbuh(2,4)=(mc*vv(2*i+n(1))+mb*vv(2*i+n(2)))*akm*abs(ya-yb)
	 Nbuh(2,5)=(mb*vv(2*i+n(2))+mc*vv(2*i+n(3)))*amn*abs(yb)	 
	 Nbuh(3,1)=(vv(2*i+n(1))-vv(2*i+n(3)))*ank*abs(ya)
	 Nbuh(3,2)=(vv(2*i+n(2))-vv(2*i+n(3)))*amn*abs(yb)
	 Nbuh(3,3)=(vv(2*i+n(2))+ma*vv(2*i+n(3)))*amn*abs(yb)+(vv(2*i+n(1))+ma*vv(2*i+n(3)))*ank*abs(ya)
	 Nbuh(3,5)=(mc*vv(2*i+n(2))+mb*vv(2*i+n(3)))*amn*abs(yb)
	 Nbuh(3,6)=(mc*vv(2*i+n(1))+mb*vv(2*i+n(3)))*ank*abs(ya)
	 
	 Nbvh(1,1)=(ma*vv(2*i+n(1))+vv(2*i+n(2)))*bkm*abs(xa-xb)+(ma*vv(2*i+n(1))+vv(2*i+n(3)))*bnk*abs(xa)
	 Nbvh(1,2)=(-vv(2*i+n(1))+vv(2*i+n(2)))*bkm*abs(xa-xb)
	 Nbvh(1,3)=(-vv(2*i+n(1))+vv(2*i+n(3)))*bnk*abs(xa)
	 Nbvh(1,4)=(mb*vv(2*i+n(1))+mc*vv(2*i+n(2)))*bkm*abs(xa-xb)
	 Nbvh(1,6)=(mb*vv(2*i+n(1))+mc*vv(2*i+n(3)))*bnk*abs(xa)
	 Nbvh(2,1)=(vv(2*i+n(1))-vv(2*i+n(2)))*bkm*abs(xa-xb)
	 Nbvh(2,2)=(vv(2*i+n(1))+ma*vv(2*i+n(2)))*bkm*abs(xa-xb)+(ma*vv(2*i+n(2))+vv(2*i+n(3)))*bmn*abs(xb)
	 Nbvh(2,3)=(-vv(2*i+n(2))+vv(2*i+n(3)))*bmn*abs(xb)
	 Nbvh(2,4)=(mc*vv(2*i+n(1))+mb*vv(2*i+n(2)))*bkm*abs(xa-xb)
	 Nbvh(2,5)=(mb*vv(2*i+n(2))+mc*vv(2*i+n(3)))*bmn*abs(xb)	 
	 Nbvh(3,1)=(vv(2*i+n(1))-vv(2*i+n(3)))*bnk*abs(xa)
	 Nbvh(3,2)=(vv(2*i+n(2))-vv(2*i+n(3)))*bmn*abs(xb)
	 Nbvh(3,3)=(vv(2*i+n(2))+ma*vv(2*i+n(3)))*bmn*abs(xb)+(vv(2*i+n(1))+ma*vv(2*i+n(3)))*bnk*abs(xa)
	 Nbvh(3,5)=(mc*vv(2*i+n(2))+mb*vv(2*i+n(3)))*bmn*abs(xb)
	 Nbvh(3,6)=(mc*vv(2*i+n(1))+mb*vv(2*i+n(3)))*bnk*abs(xa)

	 Nbu(1,1)=(ma*vv(n(1))-vv(n(2))+mb*vv(n(4)))*akm*abs(ya-yb)+(ma*vv(n(1))-vv(n(3))+mb*vv(n(6)))*ank*abs(ya)
	 Nbu(1,2)=(vv(n(1))+vv(n(2))+mc*vv(n(4)))*akm*abs(ya-yb)
	 Nbu(1,3)=(vv(n(1))+vv(n(3))+mc*vv(n(6)))*ank*abs(ya)
	 Nbu(2,1)=(vv(n(1))+vv(n(2))+mc*vv(n(4)))*akm*abs(ya-yb)
	 Nbu(2,2)=(-vv(n(1))+ma*vv(n(2))+mb*vv(n(4)))*akm*abs(ya-yb)+(ma*vv(n(2))-vv(n(3))+mb*vv(n(5)))*amn*abs(yb)
	 Nbu(2,3)=(vv(n(2))+vv(n(3))+mc*vv(n(5)))*amn*abs(yb)
	 Nbu(3,1)=(vv(n(1))+vv(n(3))+mc*vv(n(6)))*ank*abs(ya)
	 Nbu(3,2)=(vv(n(2))+vv(n(3))+mc*vv(n(5)))*amn*abs(yb)
	 Nbu(3,3)=(-vv(n(2))+ma*vv(n(3))+mb*vv(n(5)))*amn*abs(yb)+(-vv(n(1))+ma*vv(n(3))+mb*vv(n(6)))*ank*abs(ya)

	 Nbv(1,1)=(ma*vv(i+n(1))-vv(i+n(2))+mb*vv(i+n(4)))*bkm*abs(xa-xb)+(ma*vv(i+n(1))-vv(i+n(3))+mb*vv(i+n(6)))*bnk*abs(xa)
	 Nbv(1,2)=(vv(i+n(1))+vv(i+n(2))+mc*vv(i+n(4)))*bkm*abs(xa-xb)
	 Nbv(1,3)=(vv(i+n(1))+vv(i+n(3))+mc*vv(i+n(6)))*bnk*abs(xa)
	 Nbv(2,1)=(vv(i+n(1))+vv(i+n(2))+mc*vv(i+n(4)))*bkm*abs(xa-xb)
	 Nbv(2,2)=(-vv(i+n(1))+ma*vv(i+n(2))+mb*vv(i+n(4)))*bkm*abs(xa-xb)+(ma*vv(i+n(2))-vv(i+n(3))+mb*vv(i+n(5)))*bmn*abs(xb)
	 Nbv(2,3)=(vv(i+n(2))+vv(i+n(3))+mc*vv(i+n(5)))*bmn*abs(xb)
	 Nbv(3,1)=(vv(i+n(1))+vv(i+n(3))+mc*vv(i+n(6)))*bnk*abs(xa)
	 Nbv(3,2)=(vv(i+n(2))+vv(i+n(3))+mc*vv(i+n(5)))*bmn*abs(xb)
	 Nbv(3,3)=(-vv(i+n(2))+ma*vv(i+n(3))+mb*vv(i+n(5)))*bmn*abs(xb)+(-vv(i+n(1))+ma*vv(i+n(3))+mb*vv(i+n(6)))*bnk*abs(xa)

    !C�lculo de los valores de otras matrices elementales (4 puntos de integraci�n)
    !------------------------------------------------------------------------------
	do uu=1,4			   
	 Mi(1)=luno(uu)*(2.0_8*luno(uu)-1.0_8)
	 Mi(2)=ldos(uu)*(2.0_8*ldos(uu)-1.0_8)
	 Mi(3)=(1.0_8-luno(uu)-ldos(uu))*(1.0_8-2.0_8*luno(uu)-2.0_8*ldos(uu))
	 Mi(4)=4.0_8*luno(uu)*ldos(uu)
	 Mi(5)=4.0_8*ldos(uu)*(1.0_8-luno(uu)-ldos(uu))
	 Mi(6)=4.0_8*luno(uu)*(1.0_8-luno(uu)-ldos(uu))

	 Mp(1)=luno(uu)
	 Mp(2)=ldos(uu)
	 Mp(3)=1.0_8-luno(uu)-ldos(uu)

	 Mpx(1)=(y(n(2))-y(n(3)))/jac
	 Mpx(2)=(y(n(3))-y(n(1)))/jac
	 Mpx(3)=(y(n(1))-y(n(2)))/jac
	 
	 Mpy(1)=(x(n(3))-x(n(2)))/jac
	 Mpy(2)=(x(n(1))-x(n(3)))/jac
	 Mpy(3)=(x(n(2))-x(n(1)))/jac

	 Mu=dot_product(vv(n),Mi)   
	 Mv=dot_product(vv(n+i),Mi)

	 do ui=1,3					   
  	 do uj=1,3
     !Posici�n 3,3. Caja Dyx.
	 Ndyx(uj,ui)=Ndyx(uj,ui)+(Mpx(uj)*Mp(ui)*Mu+Mpy(uj)*Mp(ui)*Mv)*ja(uu)*jac 
	 enddo
	 enddo 		   	  
    enddo

	!C�lculo de los valores de otras matrices elementales (7 puntos de integraci�n)
    !------------------------------------------------------------------------------
    do uu=1,7
	 Mi(1)=lunot(uu)*(2.0_8*lunot(uu)-1.0_8)
	 Mi(2)=ldost(uu)*(2.0_8*ldost(uu)-1.0_8)
	 Mi(3)=(1.0_8-lunot(uu)-ldost(uu))*(1.0_8-2.0_8*lunot(uu)-2.0_8*ldost(uu))
	 Mi(4)=4.0_8*lunot(uu)*ldost(uu)
	 Mi(5)=4.0_8*ldost(uu)*(1.0_8-lunot(uu)-ldost(uu))
	 Mi(6)=4.0_8*lunot(uu)*(1.0_8-lunot(uu)-ldost(uu))

	 Mix(1)=((4.0_8*lunot(uu)-1.0_8)*yb)/jac
	 Mix(2)=((1.0_8-4.0_8*ldost(uu))*ya)/jac
	 Mix(3)=((4.0_8*lunot(uu)+4.0_8*ldost(uu)-3.0_8)*(yb-ya))/jac
	 Mix(4)=(4.0_8*ldost(uu)*yb-4.0_8*lunot(uu)*ya)/jac
	 Mix(5)=((8.0_8*ldost(uu)+4.0_8*lunot(uu)-4.0_8)*ya-4.0_8*ldost(uu)*yb)/jac
	 Mix(6)=((4.0_8-4.0_8*ldost(uu)-8.0_8*lunot(uu))*yb+4.0_8*lunot(uu)*ya)/jac
	
	 Miy(1)=((1.0_8-4.0_8*lunot(uu))*xb)/jac
	 Miy(2)=((4.0_8*ldost(uu)-1.0_8)*xa)/jac
	 Miy(3)=((4.0_8*lunot(uu)+4.0_8*ldost(uu)-3.0_8)*(xa-xb))/jac
	 Miy(4)=(4.0_8*lunot(uu)*xa-4.0_8*ldost(uu)*xb)/jac
	 Miy(5)=(4.0_8*ldost(uu)*xb+(4.0_8-4.0_8*lunot(uu)-8.0_8*ldost(uu))*xa)/jac
	 Miy(6)=((8.0_8*lunot(uu)+4.0_8*ldost(uu)-4.0_8)*xb-4.0_8*lunot(uu)*xa)/jac

	 Mux=dot_product(vv(n),Mix)   
	 Mvx=dot_product(vv(n+i),Mix)
	 Muy=dot_product(vv(n),Miy)   
	 Mvy=dot_product(vv(n+i),Miy)
	 
     do ui=1,6
  	 do uj=1,6 	         
	 !Posici�n 1,1. Cajas DCUA y t�rmino de fricci�n (Manning).
	 Ncua(uj,ui)=Ncua(uj,ui)+Mi(uj)*(Mi(ui)*Mux)*jat(uu)*jac					 
	 Nfx(uj,ui)=Nfx(uj,ui)+(Mi(uj)*Mi(ui)*sqrt(Mui**2+Mvi**2)*s)*jat(uu)*jac  		         
	 !Posici�n 1,2. Cajas DCVA. 
	 Ncva(uj,ui)=Ncva(uj,ui)+Mi(uj)*(Mi(ui)*Muy)*jat(uu)*jac			         
	 !Posici�n 2,1. Cajas DCUB. 
	 Ncub(uj,ui)=Ncub(uj,ui)+Mi(uj)*(Mi(ui)*Mvx)*jat(uu)*jac			         
	 !Posici�n 2,2. Cajas DCVB y t�rmino de fricci�n (Manning).
	 Ncvb(uj,ui)=Ncvb(uj,ui)+Mi(uj)*(Mi(ui)*Mvy)*jat(uu)*jac					 
	 Nfy(uj,ui)=Nfy(uj,ui)+(Mi(uj)*Mi(ui)*sqrt(Mui**2+Mvi**2)*s)*jat(uu)*jac	 	         
	 enddo
	 enddo	   	  
    enddo
	 call suma(Ncua+Nfx*9.81_8-Nbnux*nu/6.0_8,n,0,0,3*i,6,6,poi)
	 call suma(Ncva,n,0,i,3*i,6,6,podosi)	 
	 call suma(Nbhx*9.81_8,n,0,2*i,3*i,6,3,potresi)	 
	 call suma(Ncub,n,i,0,3*i,6,6,poi)
	 call suma(Ncvb+Nfy*9.81_8-Nbnux*nu/6.0_8,n,i,i,3*i,6,6,podosi)	 
	 call suma(Nbhy*9.81_8,n,i,2*i,3*i,6,3,potresi)
	 call suma(Nbuh/60.0_8,n,2*i,0,3*i,3,6,poi)
     call suma(Nbvh/60.0_8,n,2*i,i,3*i,3,6,podosi)	 
	 call suma(-Ndyx+Nbu/60.0_8+Nbv/60.0_8,n,2*i,2*i,3*i,3,3,potresi)	 
 enddo
 close(1)
!A continuaci�n se siguen los procesos que se siguen en la subrutina reducciondelsistema.

!Escritura de la matriz del sistema en formato MSR (Modified Sparse Row)
!-----------------------------------------------------------------------
!Se trata de un almacenamiento indexado por filas con dos vectores.
inc=3*i
!Generaci�n de vectores con formato MSR.
call orden (inc,k)
ndim=ndim-k

!Imposici�n de condiciones de contorno sobre el sistema
!------------------------------------------------------ 
cc=0
!S�lo se almacenan en 'vic' CC nulas en las filas y columnas a eliminar por tenerse una discretizaci�n menor.
do u=1,i
 if (v(u).eq.1) then
 vic(u+inc-i)=0.0_8
 endif
enddo 

!Trozo diferente al de la subrutina reducci�ndelsistema. Aqu� no es necesario modificar el t�rmino independiente con las CC almacenadas en 'vic' 
!como en la subrutina reducci�ndelsistema, pues se introducir�n valores nulos como valores conocidos de las inc�gnitas. 
!Por tanto, no es necesario calcular vi y simplemente se reduce el orden del t�rmino independiente (se eliminan filas).   
do u=1,inc
 !El sistema se reduce en aquellos nodos donde haya CC (valor distinto de ra�z de dos) no nula o CC nula.
 if (vic(u).eq.sqrt(2.0_8)) then 
 cc=cc+1
 vector(cc)=vector(u) 
 vb(cc)=vb(u)
 endif
enddo

!Se reduce el orden de la matriz del sistema (se eliminan filas y columnas). Se mantiene el formato MSR durante este proceso.
call reduccion(ndim,vic,inc)	  

write(6,*) ' ' 
write(6,*)'Longitud del vector matriz dispersa:',ndim

!Escritura de la matriz del sistema en formato CSC (Compressed Sparse Column)
!----------------------------------------------------------------------------
!Se trata de un almacenamiento indexado por columnas con tres vectores.
if (bcg.eq.'no') then
 !Se deja el formato MSR si se va a calcular con precondicionador diagonal (requerimientos de la programaci�n del solver).
 !Vectores 'ita,sa' con coeficientes hasta 'ndim' y con formato MSR. Aunque su dimensi�n es mucho mayor se utilizar�n estos vectores y s�lo 
 !esos coeficientes evitando generar nuevos vectores de dimensi�n 'ndim' (copiando 'ita,sa' a 'ila,la' tras allocate(ila(ndim),la(ndim))).  
 !Los vectores 'cia,ca' no son necesarios.
 deallocate(cia,ca)
else
 !Se pasa a formato CSC si se va a calcular con precondicionador LU (requerimientos de la programaci�n del solver).
 allocate(cja(cc+1))
 !Inicializaci�n de variables.
 do u=1,3*i
 ck(u)=0
 enddo
 cja(1)=1
 do u=1,cc
 cja(u+1)=0
 enddo
 !C�lculo del n�mero de coeficientes por columna.
 do u=cc+2,ndim
 cja(ita(u)+1)=cja(ita(u)+1)+1
 enddo
 !Configuraci�n final de 'cja' (considerando el coeficiente de la diagonal) y escritura de la diagonal en 'cia,ca'.
 !Adem�s se utiliza el contador 'ck'.
 do u=1,cc						
 cja(u+1)=cja(u+1)+cja(u)+1
 ca(cja(u))=sa(u)
 cia(cja(u))=u
 ck(u)=cja(u)+1
 enddo
 !Escritura del resto de coeficientes en 'cia,ca'.
 do u=1,cc										   
  do w=ita(u),ita(u+1)-1
  !En 'w' est� la columna y en 'ck(w)' la posici�n preparada para ese coeficiente.
  cia(ck(ita(w)))=u
  ca(ck(ita(w)))=sa(w)
  ck(ita(w))=ck(ita(w))+1
  enddo
 enddo
 !Vectores 'cia,ca' con coeficientes hasta 'ndim'-1 y con formato CSC. Aunque su dimensi�n es mucho mayor se utilizar�n estos vectores y s�lo esos 
 !coeficientes evitando generar nuevos vectores de dimensi�n 'ndim'-1 (copiando 'ita,sa' a 'dia,da' tras allocate(dia(ndim-1),da(ndim-1),cja(cc+1))). 
 !Los vectores 'ita,sa' no son necesarios.
 deallocate(ita,sa)
endif
deallocate(poi,podosi,potresi,ck)
end
	 
!-------------------------------------------------------------------------------------------------------------------------------------
!Subrutina CAJASASUBT
!A continuaci�n se crea la caja As que aparece en la ecuaci�n para flujo subterr�neo (con esta ecuaci�n la discretizaci�n es 
!siempre con elementos lineales, adem�s todas las cajas ir�n en la misma posici�n). 
!'Nn' son las matrices elementales de esta caja. En 'vv' est�n los niveles fre�ticos soluci�n de la iteraci�n anterior de los que 
!se obtiene el espesor del que dependen estas cajas. Por ello se trata de una caja no lineal del sistema.
!Tambi�n se da la condici�n de espesor m�nimo en caso de que el espesor obtenido sea negativo.
!------------------------------------------------------------------------------------------------------------------------------------- 
subroutine cajasasubt(i,j,x,y,zp,vv,kix,kiy,ag)
use elemental
integer*4, dimension(:),allocatable::po
integer*4 i,j
real*8, dimension(:),allocatable::kxx,kyy,kxy
real*8 x(i),y(i),zp(i),Nn(3,3),Mpx(3),Mpy(3),Mh,P,kex,key,kexy,vv(i),kix(i),kiy(i),ag(i),a,b

allocate(po(i),kxx(i),kyy(i),kxy(i))
   
39     format(4/,A80)
40     format(6X,6(X,I5))

!Inizializaci�n previa de variables:
!-----------------------------------
luno=(/0.5_8, 0.0_8, 0.5_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8/)		   
ldos=(/0.5_8, 0.5_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8/)
ja=(/1.0_8/6.0_8, 1.0_8/6.0_8, 1.0_8/6.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8/)
do u=1,i
po(u)=pos(u)
enddo
!C�lculo de las conductividades kxx, kyy y kxy a partir de las conductividades en las direcciones principales y el 
!�ngulo de anisotrop�a.
do u=1,i
kxx(u)=kix(u)*(cos(ag(u))**2)+kiy(u)*(sin(ag(u))**2)
kyy(u)=kix(u)*(sin(ag(u))**2)+kiy(u)*(cos(ag(u))**2)
kxy(u)=(kix(u)-kiy(u))*sin(ag(u))*cos(ag(u))
enddo

!C�lculo de los valores para cada elemento
!-----------------------------------------
open(unit=3,file='C:\mallasub.txt',status='old')
read(3,39)ac
 do u=1,j
     read(3,40) n(1),n(4),n(2),n(5),n(3),n(6)
  
     xa=x(n(1))-x(n(3)) 
	 xb=x(n(2))-x(n(3)) 
	 ya=y(n(1))-y(n(3))
	 yb=y(n(2))-y(n(3))
	 jac=abs(xa*yb-xb*ya)

	 do uu=1,3
     do uuu=1,3
     Nn(uu,uuu)=0.0_8
     enddo
     enddo

	 do uu=1,3
	 Mpx(1)=(y(n(2))-y(n(3)))/jac
	 Mpx(2)=(y(n(3))-y(n(1)))/jac
	 Mpx(3)=(y(n(1))-y(n(2)))/jac
	 
	 Mpy(1)=(x(n(3))-x(n(2)))/jac
	 Mpy(2)=(x(n(1))-x(n(3)))/jac
	 Mpy(3)=(x(n(2))-x(n(1)))/jac  

	 !Se calculan valores en cada punto de integraci�n del nivel fre�tico, el sustrato impermeable y conductividades.
	 Mh=luno(uu)*vv(n(1))+ldos(uu)*vv(n(2))+(1.0_8-luno(uu)-ldos(uu))*vv(n(3))
	 P=luno(uu)*zp(n(1))+ldos(uu)*zp(n(2))+(1.0_8-luno(uu)-ldos(uu))*zp(n(3))
	 kex=luno(uu)*kxx(n(1))+ldos(uu)*kxx(n(2))+(1.0_8-luno(uu)-ldos(uu))*kxx(n(3))
	 key=luno(uu)*kyy(n(1))+ldos(uu)*kyy(n(2))+(1.0_8-luno(uu)-ldos(uu))*kyy(n(3))
	 kexy=luno(uu)*kxy(n(1))+ldos(uu)*kxy(n(2))+(1.0_8-luno(uu)-ldos(uu))*kxy(n(3))
	 
	 !Condici�n de espesor m�nimo. En la ecuaci�n de agua subterr�nea al igual que las ecuaciones de aguas someras no se pueden introducir valores
	 !negativos de espesor. Aqu� no se modifica el valor de nivel fre�tico de la iteraci�n anterior, se reemplaza el espesor en el punto de 
	 !integraci�n (Mh-P) por 0.001 si se obtiene un valor negativo. Es equivalente a que se modifique la conductividad en estos elementos de forma
	 !que la transmisividad sea positiva.
	 if ((Mh-P).le.0.0) then
	 Mh=P+0.001_8
	 endif 
	 		
	  do ui=1,3
  	  do uj=1,3
	  !Posicion 1,1. Cajas Asx+Asy.
	  a=kex*(Mh-P)*Mpx(ui)+kexy*(Mh-P)*Mpy(ui)
	  b=kexy*(Mh-P)*Mpx(ui)+key*(Mh-P)*Mpy(ui)
	  Nn(uj,ui)=Nn(uj,ui)+(Mpx(uj)*a+Mpy(uj)*b)*ja(uu)*jac
	  enddo
	  enddo	   
    enddo
	call suma(Nn,n,0,0,i,3,3,po)
enddo
close(3)
deallocate(po,kxx,kyy,kxy)
end
 
!------------------------------------------------------------------------------------------------------------------------------
!Subrutina TIMESUBT
!A continuaci�n se calculan coeficientes de las matrices elementales que forman la caja Ns, t�picamente llamada matriz de masa.
!'Nn' son las matrices elementales de esta caja. Estas cajas aparecen al discretizar la ecuaci�n para flujo subterr�neo.
!S�lo entra aqu� en caso de resolver el problema no estacionario.
!------------------------------------------------------------------------------------------------------------------------------
subroutine timesubt(i,j,x,y,nd,At)
use elemental
integer*4, dimension(:),allocatable::po
integer*4 i,j				  
real*8 x(i),y(i),Nn(3,3),Mp(3),nd(i),Np,At

allocate(po(i))

39     format(4/,A80)
40     format(6X,6(X,I5))
					  
!Inizializaci�n previa de variables:
!-----------------------------------
luno=(/1.0_8/3.0_8, 0.6_8, 0.2_8, 0.2_8, 0.0_8, 0.0_8, 0.0_8/)		   
ldos=(/1.0_8/3.0_8, 0.2_8, 0.6_8, 0.2_8, 0.0_8, 0.0_8, 0.0_8/)
ja=(/-27.0_8/96.0_8, 25.0_8/96.0_8, 25.0_8/96.0_8, 25.0_8/96.0_8, 0.0_8, 0.0_8, 0.0_8/)
do u=1,i
po(u)=pos(u)
enddo

!C�lculo de los valores para cada elemento
!-----------------------------------------
open(unit=3,file='C:\mallasub.txt',status='old')
read(3,39)ac	 
 do u=1,j
     read(3,40) n(1),n(4),n(2),n(5),n(3),n(6)
  
	 xa=x(n(1))-x(n(3)) 
	 xb=x(n(2))-x(n(3)) 
	 ya=y(n(1))-y(n(3))
	 yb=y(n(2))-y(n(3))  	 
	 jac=abs(xa*yb-xb*ya)
	 	 
	 do uu=1,3
     do uuu=1,3
     Nn(uu,uuu)=0.0_8
     enddo
     enddo

	 do uu=1,4			   
	 Mp(1)=luno(uu)
	 Mp(2)=ldos(uu)
	 Mp(3)=1.0_8-luno(uu)-ldos(uu)
	 
	 !Se calculan valores en cada punto de integraci�n de la porosidad efectiva.
	 Np=Mp(1)*nd(n(1))+Mp(2)*nd(n(2))+Mp(3)*nd(n(3))

	  do ui=1,3
  	  do uj=1,3
      !Posici�n 1,1. Caja Ns.
	  Nn(uj,ui)=Nn(uj,ui)+(Mp(uj)*Np*Mp(ui))*ja(uu)*jac
	  enddo
	  enddo	  	 
	 enddo
	 call suma(Nn/At,n,0,0,i,3,3,po)
 enddo
 close(3)
 deallocate(po)
end

!--------------------------------------------------------------------------------------------------------------------------
!Subrutina TIMESUBTCONC
!Es una alternativa a la subrutina timesubt. Se diferencia de la anterior en que se concentra la masa en la diagonal.	
!En otras palabras, se suman todos los t�rminos de cada fila en el t�rmino de la diagonal.
!--------------------------------------------------------------------------------------------------------------------------
subroutine timesubtconc(i,j,x,y,nd,At)
use elemental
integer*4, dimension(:),allocatable::po
integer*4 i,j				  
real*8 x(i),y(i),Nn(3,3),Mp(3),nd(i),Np,At

allocate(po(i))

39     format(4/,A80)
40     format(6X,6(X,I5))
				
!Inizializaci�n previa de variables:
!-----------------------------------
luno=(/0.5_8, 0.0_8, 0.5_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8/)		   
ldos=(/0.5_8, 0.5_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8/)
ja=(/1.0_8/6.0_8, 1.0_8/6.0_8, 1.0_8/6.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8/)
do u=1,i
po(u)=pos(u)
enddo

!C�lculo de los valores para cada elemento
!-----------------------------------------
open(unit=3,file='C:\mallasub.txt',status='old')
read(3,39)ac	 
 do u=1,j
     read(3,40) n(1),n(4),n(2),n(5),n(3),n(6)
  
	 xa=x(n(1))-x(n(3)) 
	 xb=x(n(2))-x(n(3)) 
	 ya=y(n(1))-y(n(3))
	 yb=y(n(2))-y(n(3))  	 
	 jac=abs(xa*yb-xb*ya)
	 	 
	 do uu=1,3
     do uuu=1,3
     Nn(uu,uuu)=0.0_8
     enddo
     enddo

	 do uu=1,3
	 Mp(1)=luno(uu)
	 Mp(2)=ldos(uu)
	 Mp(3)=1.0_8-luno(uu)-ldos(uu)
	 
	 Np=Mp(1)*nd(n(1))+Mp(2)*nd(n(2))+Mp(3)*nd(n(3))
 
  	  do uj=1,3
      !Posici�n 1,1. Caja Ns concentrada.
	  !La suma de las funciones de forma en cualquier punto del elemento ser� igual a la unidad. 
	  Nn(uj,uj)=Nn(uj,uj)+(Mp(uj)*Np)*ja(uu)*jac
	  enddo	  	 
	 enddo
	 call suma(Nn/At,n,0,0,i,3,3,po)
 enddo
 close(3)
 deallocate(po)
end

!---------------------------------------------------------------------------------------------------------------------------------------------
!Subrutinas para el c�lculo de coeficientes de las integrales en el dominio (formulaci�n Petrov-Galerkin).
!M�todo de estabilizaci�n SUPG/PSPG con grad-div para las ecuaciones para flujo superficial. Par�metros de estabilizaci�n fueron desarrollados 
!por Gelhard para elementos P2-P1 para mallas uniformes considerando las ecuaciones de Navier-Stokes 2D en funci�n de la presi�n cinem�tica y 
!sin considerar las integrales de contorno de las tensiones viscosas. La estabilizaci�n no est� programada para el m�todo de Newton.  
!Para elementos cuadr�ticos 'Mixx' ser�n las derivadas segundas en x de las funciones de forma y 'Miyy' las derivadas segundas en y de las 
!funciones de forma.
!---------------------------------------------------------------------------------------------------------------------------------------------
!Subrutina TIMEASUPGNS
!Para las ecuaciones de Navier-Stokes 2D. No se aplica como peso el t�rmino convectivo (despreciable) ni el t�rmino temporal.
!Se utiliza como longitud elemental el di�metro de la circunferencia que inscribe al elemento y unos par�metros similares a los de 
!Gelhard (aqu� las ecuaciones est�n escritas en funci�n del calado).
!La estabilizaci�n contiene en las dos primeras ecuaciones una ec. din�mica + ec. din�mica (con t�rmino supg)+ ec. cont (con t�rmino lsic)
!y como �ltima ecuaci�n la ec. continuidad + ec. din�mica (con peso pspg).
!Aqu� se calculan las matrices elementales de los t�rminos de masa estabilizados. Se utilizar�an al resolver de forma transitoria (en el 
!esquema semi-impl�cito siempre) las ecuaciones, pero de momento no se utilizan.
!---------------------------------------------------------------------------------------------------------------------------------------------
subroutine timeasupgns (i,j,x,y,vv,At)
use elemental
integer*4, dimension(:),allocatable::poi,podosi
integer*4 i,j
real*8 x(i),y(i),Mi(6),Mix(6),Miy(6),Mpx(3),Mpy(3),Nt(6,6),Ntpx(3,6),Ntpy(3,6)	 
real*8 At,lsic,Mu,Mv,lunoq(13),ldosq(13),jaq(13),vv(3*i),cx,ccx,h,Muu,Mvv,a,b,c 

allocate(poi(3*i),podosi(3*i))

39     format(4/,A80)
40     format(6X,6(X,I5))
  
!Inizializaci�n previa de variables:
!-----------------------------------
lunoq=(/0.0651301029002_8, 0.8697397941956_8, 0.0651301029002_8, 0.3128654960049_8, 0.6384441885698_8, 0.0486903154253_8, 	 &
0.6384441885698_8, 0.3128654960049_8, 0.0486903154253_8, 0.2603459660790_8, 0.4793080678419_8, 0.2603459660790_8, 1.0_8/3.0_8/)		   
ldosq=(/0.0651301029002_8, 0.0651301029002_8, 0.8697397941956_8, 0.0486903154253_8, 0.3128654960049_8, 0.6384441885698_8, 	 &
0.0486903154253_8, 0.6384441885698_8, 0.3128654960049_8, 0.2603459660790_8, 0.2603459660790_8, 0.4793080678419_8, 1.0_8/3.0_8/)
jaq=(/0.0266736178004_8, 0.0266736178004_8, 0.0266736178004_8, 0.03855688044515_8, 0.03855688044515_8, 0.03855688044515_8, 	 &
0.03855688044515_8, 0.03855688044515_8, 0.03855688044515_8, 0.0878076287166_8, 0.0878076287166_8, 0.0878076287166_8, 		 &
-0.0747850222335_8/)
do u=1,3*i
poi(u)=posi(u)
podosi(u)=posdosi(u)
enddo
lsic=1.0_8
!Introduciendo a=9.81 se considerar�a la gravedad en el peso PSPG tal y como deber�a de ser si se utiliza el operador de la ecuaci�n
!diferencial como peso. Introduciendo a=0.0 se considerar�a la estabilizaci�n SUPG con grad-div. En cualquiera de los dos casos se 
!obtendr�an peores soluciones.
!Se puede entender que se aplica la gravedad en el peso y un tercer par�metro ca=(h**2.0_8)/9.81_8 en vez de ccx donde se use 'a'
!de forma que los par�metros para SUPG y PSPG son diferentes.
a=1.0_8

!C�lculo de los valores para cada elemento
!-----------------------------------------
open(unit=1,file='C:\malla.txt',status='old')
read(1,39)ac
  do u=1,j
     read(1,40) n(1),n(4),n(2),n(5),n(3),n(6)

	 xa=x(n(1))-x(n(3)) 
	 xb=x(n(2))-x(n(3)) 
	 ya=y(n(1))-y(n(3))
	 yb=y(n(2))-y(n(3))  
	 jac=abs(xa*yb-xb*ya)
	 	 
	 do uu=1,6
     do uuu=1,6
	 Nt(uu,uuu)=0.0_8
     enddo
     enddo
	 do uu=1,3
     do uuu=1,6
	 Ntpx(uu,uuu)=0.0_8
	 Ntpy(uu,uuu)=0.0_8
     enddo
     enddo
     
	 Mpx(1)=(y(n(2))-y(n(3)))/jac
	 Mpx(2)=(y(n(3))-y(n(1)))/jac
	 Mpx(3)=(y(n(1))-y(n(2)))/jac
	 
	 Mpy(1)=(x(n(3))-x(n(2)))/jac
	 Mpy(2)=(x(n(1))-x(n(3)))/jac
	 Mpy(3)=(x(n(2))-x(n(1)))/jac
	 
	 !Valores medios de velocidad en el elemento (valores en el centro del elemento).
	 Muu=-(vv(n(1))+vv(n(2))+vv(n(3)))/9.0_8+(vv(n(4))+vv(n(5))+vv(n(6)))*4.0_8/9.0_8   
	 Mvv=-(vv(n(1)+i)+vv(n(2)+i)+vv(n(3)+i))/9.0_8+(vv(n(4)+i)+vv(n(5)+i)+vv(n(6)+i))*4.0_8/9.0_8

	 !C�lculo del diametro de la circunferencia inscrita:
	 a=sqrt((x(n(1))-x(n(2)))**2.0_8+(y(n(1))-y(n(2)))**2.0_8)
	 b=sqrt((x(n(2))-x(n(3)))**2.0_8+(y(n(2))-y(n(3)))**2.0_8)
	 c=sqrt((x(n(3))-x(n(1)))**2.0_8+(y(n(3))-y(n(1)))**2.0_8)
	 h=2.0_8*jac/(a+b+c)

	 !Par�metros de estabilizaci�n (del tipo de los de Tobias Geldar) para elementos que cumplen la condici�n LBB. Se podr�an escribir 
	 !como cx=9.81_8, ccx=(h**2.0_8). Mismos par�metros para SUPG y PSPG (en ccx).
	 cx=9.81_8*lsic
	 ccx=9.81_8*(h**2.0_8)/cx

	 !Integraci�n gauss de quinto orden con trece puntos de integraci�n. Polinomios de grado superior a grado 7 quedar�an 
	 !integrados con error.
	 do uu=1,13
	 
	 Mi(1)=lunoq(uu)*(2.0_8*lunoq(uu)-1.0_8)
	 Mi(2)=ldosq(uu)*(2.0_8*ldosq(uu)-1.0_8)
	 Mi(3)=(1.0_8-lunoq(uu)-ldosq(uu))*(1.0_8-2.0_8*lunoq(uu)-2.0_8*ldosq(uu))
	 Mi(4)=4.0_8*lunoq(uu)*ldosq(uu)
	 Mi(5)=4.0_8*ldosq(uu)*(1.0_8-lunoq(uu)-ldosq(uu))
	 Mi(6)=4.0_8*lunoq(uu)*(1.0_8-lunoq(uu)-ldosq(uu))

	 Mix(1)=((4.0_8*lunoq(uu)-1.0_8)*yb)/jac
	 Mix(2)=((1.0_8-4.0_8*ldosq(uu))*ya)/jac
	 Mix(3)=((4.0_8*lunoq(uu)+4.0_8*ldosq(uu)-3.0_8)*(yb-ya))/jac
	 Mix(4)=(4.0_8*ldosq(uu)*yb-4.0_8*lunoq(uu)*ya)/jac
	 Mix(5)=((8.0_8*ldosq(uu)+4.0_8*lunoq(uu)-4.0_8)*ya-4.0_8*ldosq(uu)*yb)/jac
	 Mix(6)=((4.0_8-4.0_8*ldosq(uu)-8.0_8*lunoq(uu))*yb+4.0_8*lunoq(uu)*ya)/jac
	
	 Miy(1)=((1.0_8-4.0_8*lunoq(uu))*xb)/jac
	 Miy(2)=((4.0_8*ldosq(uu)-1.0_8)*xa)/jac
	 Miy(3)=((4.0_8*lunoq(uu)+4.0_8*ldosq(uu)-3.0_8)*(xa-xb))/jac
	 Miy(4)=(4.0_8*lunoq(uu)*xa-4.0_8*ldosq(uu)*xb)/jac
	 Miy(5)=(4.0_8*ldosq(uu)*xb+(4.0_8-4.0_8*lunoq(uu)-8.0_8*ldosq(uu))*xa)/jac
	 Miy(6)=((8.0_8*lunoq(uu)+4.0_8*ldosq(uu)-4.0_8)*xb-4.0_8*lunoq(uu)*xa)/jac
 
	 Mu=dot_product(vv(n),Mi)   
	 Mv=dot_product(vv(n+i),Mi)

	  do ui=1,6
  	  do uj=1,6 
	  !Posici�n 1,1 y 2,2. Cajas de las matrices de masa M afectadas por el peso SUPG. 
	  Nt(uj,ui)=Nt(uj,ui)+ccx*(Mu*Mix(uj)+Mv*Miy(uj))*Mi(ui)*jaq(uu)*jac
	  enddo								  
	  enddo							 
	  
	  do ui=1,6	
	  do uj=1,3
	  !Posicion 3,3. Cajas de las matrices de masa M afectadas por el peso PSPG.	 	  	 	  
	  Ntpx(uj,ui)=Ntpx(uj,ui)+ccx*Mpx(uj)*Mi(ui)*jaq(uu)*jac
	  Ntpy(uj,ui)=Ntpy(uj,ui)+ccx*Mpy(uj)*Mi(ui)*jaq(uu)*jac 
	  enddo
	  enddo	    
	 enddo	 
	 call suma(Nt/At,n,0,0,3*i,6,6,poi)
	 call suma(Nt/At,n,i,i,3*i,6,6,podosi)	 
	 !Sin las siguientes matrices s�lo se tendr�a la estabilizaci�n SUPG con grad-div.
	 call suma(a*Ntpx/At,n,2*i,0,3*i,3,6,poi)
	 call suma(a*Ntpy/At,n,2*i,i,3*i,3,6,podosi)
 enddo
 close(1)
 deallocate(poi,podosi)
end

!-----------------------------------------------------------------------------------------------------------------------
!Subrutina CAJASUPGNS
!Para las ecuaciones de Navier-Stokes 2D. Aqu� se calculan las matrices elementales del resto de t�rminos estabilizados.
!-----------------------------------------------------------------------------------------------------------------------
subroutine cajasupgns(i,j,x,y,vv,nu,del)
use elemental	  
integer*4, dimension(:),allocatable::poi,podosi,potresi
integer*4 i,j
real*8 x(i),y(i),Nn(6,6),Nm(6,6),Nk(6,3),Np(6,3),Mi(6),Mix(6),Miy(6),Mixx(6),Miyy(6),Mp(3),Mpx(3),Mpy(3)
real*8 Nh(3,6),Nf(3,6),Nl(3,3),Nj(3,6),Nw(3,6),Nxx(6,6),Nxy(6,6),Nyx(6,6),Nyy(6,6),lsic
real*8 Mu,Mv,lunoq(13),ldosq(13),jaq(13),vv(3*i),cx,ccx,nu,h,Muu,Mvv,a,b,c,del 

allocate(poi(3*i),podosi(3*i),potresi(3*i))

39     format(4/,A80)
40     format(6X,6(X,I5))

!Inizializaci�n previa de variables:
!-----------------------------------
lunoq=(/0.0651301029002_8, 0.8697397941956_8, 0.0651301029002_8, 0.3128654960049_8, 0.6384441885698_8, 0.0486903154253_8, 	 &
0.6384441885698_8, 0.3128654960049_8, 0.0486903154253_8, 0.2603459660790_8, 0.4793080678419_8, 0.2603459660790_8, 1.0_8/3.0_8/)		   
ldosq=(/0.0651301029002_8, 0.0651301029002_8, 0.8697397941956_8, 0.0486903154253_8, 0.3128654960049_8, 0.6384441885698_8, 	 &
0.0486903154253_8, 0.6384441885698_8, 0.3128654960049_8, 0.2603459660790_8, 0.2603459660790_8, 0.4793080678419_8, 1.0_8/3.0_8/)
jaq=(/0.0266736178004_8, 0.0266736178004_8, 0.0266736178004_8, 0.03855688044515_8, 0.03855688044515_8, 0.03855688044515_8, 	 &
0.03855688044515_8, 0.03855688044515_8, 0.03855688044515_8, 0.0878076287166_8, 0.0878076287166_8, 0.0878076287166_8, 		 &
-0.0747850222335_8/)
do u=1,3*i
poi(u)=posi(u)
podosi(u)=posdosi(u)
potresi(u)=postresi(u)
enddo
lsic=1.0_8
a=1.0_8
  
!C�lculo de los valores para cada elemento
!-----------------------------------------
open(unit=1,file='C:\malla.txt',status='old')
read(1,39)ac
  do u=1,j
     read(1,40) n(1),n(4),n(2),n(5),n(3),n(6)

	 xa=x(n(1))-x(n(3)) 
	 xb=x(n(2))-x(n(3)) 
	 ya=y(n(1))-y(n(3))
	 yb=y(n(2))-y(n(3))  
	 jac=abs(xa*yb-xb*ya)
	 	 
	 do uu=1,6
     do uuu=1,6
     Nn(uu,uuu)=0.0_8
	 Nm(uu,uuu)=0.0_8	 
	 Nxx(uu,uuu)=0.0_8
	 Nxy(uu,uuu)=0.0_8
	 Nyx(uu,uuu)=0.0_8
	 Nyy(uu,uuu)=0.0_8
     enddo
     do uuu=1,3
	 Nk(uu,uuu)=0.0_8
	 Np(uu,uuu)=0.0_8	 
     enddo
     enddo
	 do uu=1,3
     do uuu=1,6
	 Nh(uu,uuu)=0.0_8
	 Nf(uu,uuu)=0.0_8	 
     enddo
     enddo
     
	 !Se calculan antes estos valores ya que son contantes para cada elemento.
	 Mixx(1)=4.0_8*(yb**2.0_8)/(jac**2.0_8)
	 Mixx(2)=4.0_8*(ya**2.0_8)/(jac**2.0_8)
	 Mixx(3)=4.0_8*((yb-ya)**2.0_8)/(jac**2.0_8)
	 Mixx(4)=-8.0_8*ya*yb/(jac**2.0_8)
	 Mixx(5)=8.0_8*(ya*yb-ya**2.0_8)/(jac**2.0_8)
	 Mixx(6)=8.0_8*(ya*yb-yb**2.0_8)/(jac**2.0_8)

	 Miyy(1)=4.0_8*(xb**2.0_8)/(jac**2.0_8)
	 Miyy(2)=4.0_8*(xa**2.0_8)/(jac**2.0_8)
	 Miyy(3)=4.0_8*((xa-xb)**2.0_8)/(jac**2.0_8)
	 Miyy(4)=-8.0_8*xa*xb/(jac**2.0_8)
	 Miyy(5)=8.0_8*(xa*xb-xa**2.0_8)/(jac**2.0_8)
	 Miyy(6)=8.0_8*(xa*xb-xb**2.0_8)/(jac**2.0_8)

	 Mpx(1)=(y(n(2))-y(n(3)))/jac
	 Mpx(2)=(y(n(3))-y(n(1)))/jac
	 Mpx(3)=(y(n(1))-y(n(2)))/jac
	 
	 Mpy(1)=(x(n(3))-x(n(2)))/jac
	 Mpy(2)=(x(n(1))-x(n(3)))/jac
	 Mpy(3)=(x(n(2))-x(n(1)))/jac
	 
	 Muu=-(vv(n(1))+vv(n(2))+vv(n(3)))/9.0_8+(vv(n(4))+vv(n(5))+vv(n(6)))*4.0_8/9.0_8   
	 Mvv=-(vv(n(1)+i)+vv(n(2)+i)+vv(n(3)+i))/9.0_8+(vv(n(4)+i)+vv(n(5)+i)+vv(n(6)+i))*4.0_8/9.0_8
	 
	 a=sqrt((x(n(1))-x(n(2)))**2.0_8+(y(n(1))-y(n(2)))**2.0_8)
	 b=sqrt((x(n(2))-x(n(3)))**2.0_8+(y(n(2))-y(n(3)))**2.0_8)
	 c=sqrt((x(n(3))-x(n(1)))**2.0_8+(y(n(3))-y(n(1)))**2.0_8)
	 h=2.0_8*jac/(a+b+c)

	 cx=9.81_8*lsic
	 ccx=9.81_8*(h**2.0_8)/cx
	 
	 do uu=1,13
     
	 Mi(1)=lunoq(uu)*(2.0_8*lunoq(uu)-1.0_8)
	 Mi(2)=ldosq(uu)*(2.0_8*ldosq(uu)-1.0_8)
	 Mi(3)=(1.0_8-lunoq(uu)-ldosq(uu))*(1.0_8-2.0_8*lunoq(uu)-2.0_8*ldosq(uu))
	 Mi(4)=4.0_8*lunoq(uu)*ldosq(uu)
	 Mi(5)=4.0_8*ldosq(uu)*(1.0_8-lunoq(uu)-ldosq(uu))
	 Mi(6)=4.0_8*lunoq(uu)*(1.0_8-lunoq(uu)-ldosq(uu))

	 Mix(1)=((4.0_8*lunoq(uu)-1.0_8)*yb)/jac
	 Mix(2)=((1.0_8-4.0_8*ldosq(uu))*ya)/jac
	 Mix(3)=((4.0_8*lunoq(uu)+4.0_8*ldosq(uu)-3.0_8)*(yb-ya))/jac
	 Mix(4)=(4.0_8*ldosq(uu)*yb-4.0_8*lunoq(uu)*ya)/jac
	 Mix(5)=((8.0_8*ldosq(uu)+4.0_8*lunoq(uu)-4.0_8)*ya-4.0_8*ldosq(uu)*yb)/jac
	 Mix(6)=((4.0_8-4.0_8*ldosq(uu)-8.0_8*lunoq(uu))*yb+4.0_8*lunoq(uu)*ya)/jac
	
	 Miy(1)=((1.0_8-4.0_8*lunoq(uu))*xb)/jac
	 Miy(2)=((4.0_8*ldosq(uu)-1.0_8)*xa)/jac
	 Miy(3)=((4.0_8*lunoq(uu)+4.0_8*ldosq(uu)-3.0_8)*(xa-xb))/jac
	 Miy(4)=(4.0_8*lunoq(uu)*xa-4.0_8*ldosq(uu)*xb)/jac
	 Miy(5)=(4.0_8*ldosq(uu)*xb+(4.0_8-4.0_8*lunoq(uu)-8.0_8*ldosq(uu))*xa)/jac
	 Miy(6)=((8.0_8*lunoq(uu)+4.0_8*ldosq(uu)-4.0_8)*xb-4.0_8*lunoq(uu)*xa)/jac
 
	 Mp(1)=lunoq(uu)
	 Mp(2)=ldosq(uu)
	 Mp(3)=1.0_8-lunoq(uu)-ldosq(uu)

	 Mu=dot_product(vv(n),Mi)   
	 Mv=dot_product(vv(n+i),Mi)

	  do ui=1,6
  	  do uj=1,6 
	  !Posici�n 1,1 y 2,2. 
	  !Cajas no lineales C con peso SUPG. 
	  Nn(uj,ui)=Nn(uj,ui)+ccx*(Mu*Mix(uj)+Mv*Miy(uj))*(Mix(ui)*Mu+Miy(ui)*Mv)*jaq(uu)*jac	  
	  !Caja A con peso supg (sin forma d�bil)
	  Nm(uj,ui)=Nm(uj,ui)+ccx*(Mu*Mix(uj)+Mv*Miy(uj))*(Miyy(ui)+Mixx(ui))*jaq(uu)*jac	  
	  !Cajas correspondientes a la ecuaci�n de continuidad estabilizada con pesos grad-div (se usa el par�metro 'cx').
	  Nxx(uj,ui)=Nxx(uj,ui)+cx*Mix(uj)*Mix(ui)*jaq(uu)*jac
	  Nxy(uj,ui)=Nxy(uj,ui)+cx*Mix(uj)*Miy(ui)*jaq(uu)*jac	  
	  Nyx(uj,ui)=Nyx(uj,ui)+cx*Miy(uj)*Mix(ui)*jaq(uu)*jac
	  Nyy(uj,ui)=Nyy(uj,ui)+cx*Miy(uj)*Miy(ui)*jaq(uu)*jac	  
	  enddo								  
	  enddo	
	  
	  do uj=1,6
	  do ui=1,3
	  !Posicion 1,3. Caja Bx con peso SUPG.
	  !Sin forma d�bil	 	  
	  Nk(uj,ui)=Nk(uj,ui)+ccx*(Mu*Mix(uj)+Mv*Miy(uj))*Mpx(ui)*jaq(uu)*jac	  
	  !Posicion 2,3. Caja By con peso SUPG.
	  !Sin forma d�bil
	  Np(uj,ui)=Np(uj,ui)+ccx*(Mu*Mix(uj)+Mv*Miy(uj))*Mpy(ui)*jaq(uu)*jac	 
	  enddo
	  enddo
	  
	  do ui=1,6	
	  do uj=1,3
	  !Posicion 3,1. Caja C con peso SUPG.	 	  
	  Nh(uj,ui)=Nh(uj,ui)+ccx*Mpx(uj)*(Mix(ui)*Mu+Miy(ui)*Mv)*jaq(uu)*jac
	  !Posicion 3,2. Caja C con peso SUPG.	 	  
	  Nf(uj,ui)=Nf(uj,ui)+ccx*Mpy(uj)*(Mix(ui)*Mu+Miy(ui)*Mv)*jaq(uu)*jac	  	  
	  enddo
	  enddo	     
	 enddo
	 
	 !Las siguientes cajas (con pesos PSPG) no se integran ya que contienen coeficientes contantes para cada elemento.	 
	 do ui=1,6	
	  do uj=1,3
	  !Posicion 3,1. Caja A con peso PSPG.	 	  
	  Nj(uj,ui)=ccx*Mpx(uj)*(Miyy(ui)+Mixx(ui))*jac/2.0_8
	  !Posicion 3,2. Caja A con peso PSPG.	 	  
	  Nw(uj,ui)=ccx*Mpy(uj)*(Miyy(ui)+Mixx(ui))*jac/2.0_8	  	  
	  enddo
	 enddo
	 do ui=1,3	
	  do uj=1,3
	  !Posicion 3,3. Caja Bx+By con peso PSPG.	 	  
	  Nl(uj,ui)=ccx*(Mpx(uj)*Mpx(ui)+Mpy(uj)*Mpy(ui))*jac/2.0_8	 	   
	  enddo
	 enddo
	 call suma((Nn-nu*Nm+Nxx)/del,n,0,0,3*i,6,6,poi)
	 call suma(Nxy/del,n,0,i,3*i,6,6,podosi)
	 call suma(Nyx/del,n,i,0,3*i,6,6,poi)
	 call suma((Nn-nu*Nm+Nyy)/del,n,i,i,3*i,6,6,podosi)
	 call suma(9.81_8*Nk/del,n,0,2*i,3*i,6,3,potresi)
	 call suma(9.81_8*Np/del,n,i,2*i,3*i,6,3,potresi)	 
	 !Sin las siguientes matrices s�lo se tendr�a la estabilizaci�n SUPG con grad-div.
	 call suma((a*Nh-a*nu*Nj)/del,n,2*i,0,3*i,3,6,poi)
	 call suma((a*Nf-a*nu*Nw)/del,n,2*i,i,3*i,3,6,podosi)
	 call suma(a*9.81_8*Nl/del,n,2*i,2*i,3*i,3,3,potresi)
 enddo
 close(1)
 deallocate(poi,podosi,potresi)
end

!------------------------------------------------------------------------------------------------------------------------------------------  
!Subrutina TIMEASUPGAS
!Para las ecuaciones de aguas someras. No se aplica como peso el t�rmino convectivo (despreciable) ni el t�rmino de maning ni el 
!t�rmino temporal. De momento la estabilizaci�n no est� aplicada al t�rmino de lluvia (ir�a en t�rmino independiente de ec din�micas 
!con peso cx).
!Se utiliza como longitud elemental el di�metro de la circunferencia que inscribe al elemento y unos par�metros similares a los de la 
!estabilizaci�n de las ecuaciones de Navier-Stokes 2D.
!La estabilizaci�n contiene en las dos primeras ecuaciones una ec. din�mica + ec. din�mica (con t�rmino supg)+ ec. cont (con t�rmino lsic)
!y como �ltima ecuaci�n la ec. continuidad + ec. din�mica (con peso pspg).
!Aqu� se calculan las matrices elementales de los t�rminos de masa estabilizados. Se utilizar�an al resolver de forma transitoria (en el 
!esquema semi-impl�cito siempre) las ecuaciones, pero de momento no se utilizan.
!-----------------------------------------------------------------------------------------------------------------------------------------
subroutine timeasupgas (i,j,x,y,vv,At)
use elemental
integer*4, dimension(:),allocatable::poi,podosi,potresi
integer*4 i,j
real*8 x(i),y(i),Mi(6),Mix(6),Miy(6),Mp(3),Mpx(3),Mpy(3),Nt(6,6),Ntpx(3,6),Ntpy(3,6)
real*8 At,lsic,Mu,Mv,lunoq(13),ldosq(13),jaq(13),vv(3*i),cx,ccx,h,Muu,Mvv,Mhu,a,b,c,hxy,Mhx,Mhy,Ntx(6,3),Nty(6,3)

allocate(poi(3*i),podosi(3*i),potresi(3*i))

39     format(4/,A80)
40     format(6X,6(X,I5))

!Inizializaci�n previa de variables:
!-----------------------------------
lunoq=(/0.0651301029002_8, 0.8697397941956_8, 0.0651301029002_8, 0.3128654960049_8, 0.6384441885698_8, 0.0486903154253_8, 	 &
0.6384441885698_8, 0.3128654960049_8, 0.0486903154253_8, 0.2603459660790_8, 0.4793080678419_8, 0.2603459660790_8, 1.0_8/3.0_8/)		   
ldosq=(/0.0651301029002_8, 0.0651301029002_8, 0.8697397941956_8, 0.0486903154253_8, 0.3128654960049_8, 0.6384441885698_8, 	 &
0.0486903154253_8, 0.6384441885698_8, 0.3128654960049_8, 0.2603459660790_8, 0.2603459660790_8, 0.4793080678419_8, 1.0_8/3.0_8/)
jaq=(/0.0266736178004_8, 0.0266736178004_8, 0.0266736178004_8, 0.03855688044515_8, 0.03855688044515_8, 0.03855688044515_8, 	 &
0.03855688044515_8, 0.03855688044515_8, 0.03855688044515_8, 0.0878076287166_8, 0.0878076287166_8, 0.0878076287166_8, 		 &
-0.0747850222335_8/)
do u=1,3*i
poi(u)=posi(u)
podosi(u)=posdosi(u)
potresi(u)=postresi(u)
enddo
!Las consideraciones sobre la variable 'a' son id�nticas a las explicadas en la estabilizaci�n de la ecuaci�n de Navier-Stokes 2D.
lsic=1.0_8
a=1.0_8
  
!C�lculo de los valores para cada elemento
!-----------------------------------------
open(unit=1,file='C:\malla.txt',status='old')
read(1,39)ac
  do u=1,j
     read(1,40) n(1),n(4),n(2),n(5),n(3),n(6)

	 xa=x(n(1))-x(n(3)) 
	 xb=x(n(2))-x(n(3)) 
	 ya=y(n(1))-y(n(3))
	 yb=y(n(2))-y(n(3))  
	 jac=abs(xa*yb-xb*ya)
	 	 
	 do uu=1,6
     do uuu=1,6
	 Nt(uu,uuu)=0.0_8
     enddo
	 do uuu=1,3
	 Ntx(uu,uuu)=0.0_8
	 Nty(uu,uuu)=0.0_8
     enddo
     enddo
	 do uu=1,3
     do uuu=1,6
	 Ntpx(uu,uuu)=0.0_8
	 Ntpy(uu,uuu)=0.0_8
     enddo
     enddo
     
	 Mpx(1)=(y(n(2))-y(n(3)))/jac
	 Mpx(2)=(y(n(3))-y(n(1)))/jac
	 Mpx(3)=(y(n(1))-y(n(2)))/jac
	 
	 Mpy(1)=(x(n(3))-x(n(2)))/jac
	 Mpy(2)=(x(n(1))-x(n(3)))/jac
	 Mpy(3)=(x(n(2))-x(n(1)))/jac

	 !'Mhx, Mhy' son las derivadas del calado en direcci�n x y en direcci�n y evaluadas en cualquier punto de integraci�n del elemento, 
	 !a partir de los valores de calado de la iteraci�n anterior. Su valor es constante en todos los puntos.
	 Mhx=vv(2*i+n(1))*Mpx(1)+vv(2*i+n(2))*Mpx(2)+vv(2*i+n(3))*Mpx(3)
	 Mhy=vv(2*i+n(1))*Mpy(1)+vv(2*i+n(2))*Mpy(2)+vv(2*i+n(3))*Mpy(3)

	 !Valores medios de velocidad y calado en el elemento (valores en el centro del elemento).
	 Muu=-(vv(n(1))+vv(n(2))+vv(n(3)))/9.0_8+(vv(n(4))+vv(n(5))+vv(n(6)))*4.0_8/9.0_8   
	 Mvv=-(vv(n(1)+i)+vv(n(2)+i)+vv(n(3)+i))/9.0_8+(vv(n(4)+i)+vv(n(5)+i)+vv(n(6)+i))*4.0_8/9.0_8
	 Mhu=(vv(n(1)+2*i)+vv(n(2)+2*i)+vv(n(3)+2*i))/3.0_8

	 !C�lculo del diametro de la circunferencia inscrita:
	 a=sqrt((x(n(1))-x(n(2)))**2.0_8+(y(n(1))-y(n(2)))**2.0_8)
	 b=sqrt((x(n(2))-x(n(3)))**2.0_8+(y(n(2))-y(n(3)))**2.0_8)
	 c=sqrt((x(n(3))-x(n(1)))**2.0_8+(y(n(3))-y(n(1)))**2.0_8)
	 h=2.0_8*jac/(a+b+c)

	 !Par�metros de estabilizaci�n en funci�n de 'Mhu' para elementos que cumplen la condici�n LBB.	Se podr�an escribir 
	 !como cx=9.81_8/(Mhu**2.0_8), ccx=(h**2.0_8). Se utilizan los mismos par�metros para SUPG y PSPG (en ccx).  
	 cx=9.81_8*lsic/(Mhu**2.0_8)
	 ccx=9.81_8*(h**2.0_8)/(cx*(Mhu**2.0_8))
	 
	 do uu=1,13
     
	 Mi(1)=lunoq(uu)*(2.0_8*lunoq(uu)-1.0_8)
	 Mi(2)=ldosq(uu)*(2.0_8*ldosq(uu)-1.0_8)
	 Mi(3)=(1.0_8-lunoq(uu)-ldosq(uu))*(1.0_8-2.0_8*lunoq(uu)-2.0_8*ldosq(uu))
	 Mi(4)=4.0_8*lunoq(uu)*ldosq(uu)
	 Mi(5)=4.0_8*ldosq(uu)*(1.0_8-lunoq(uu)-ldosq(uu))
	 Mi(6)=4.0_8*lunoq(uu)*(1.0_8-lunoq(uu)-ldosq(uu))

	 Mix(1)=((4.0_8*lunoq(uu)-1.0_8)*yb)/jac
	 Mix(2)=((1.0_8-4.0_8*ldosq(uu))*ya)/jac
	 Mix(3)=((4.0_8*lunoq(uu)+4.0_8*ldosq(uu)-3.0_8)*(yb-ya))/jac
	 Mix(4)=(4.0_8*ldosq(uu)*yb-4.0_8*lunoq(uu)*ya)/jac
	 Mix(5)=((8.0_8*ldosq(uu)+4.0_8*lunoq(uu)-4.0_8)*ya-4.0_8*ldosq(uu)*yb)/jac
	 Mix(6)=((4.0_8-4.0_8*ldosq(uu)-8.0_8*lunoq(uu))*yb+4.0_8*lunoq(uu)*ya)/jac
	
	 Miy(1)=((1.0_8-4.0_8*lunoq(uu))*xb)/jac
	 Miy(2)=((4.0_8*ldosq(uu)-1.0_8)*xa)/jac
	 Miy(3)=((4.0_8*lunoq(uu)+4.0_8*ldosq(uu)-3.0_8)*(xa-xb))/jac
	 Miy(4)=(4.0_8*lunoq(uu)*xa-4.0_8*ldosq(uu)*xb)/jac
	 Miy(5)=(4.0_8*ldosq(uu)*xb+(4.0_8-4.0_8*lunoq(uu)-8.0_8*ldosq(uu))*xa)/jac
	 Miy(6)=((8.0_8*lunoq(uu)+4.0_8*ldosq(uu)-4.0_8)*xb-4.0_8*lunoq(uu)*xa)/jac
 
	 Mp(1)=lunoq(uu)
	 Mp(2)=ldosq(uu)
	 Mp(3)=1.0_8-lunoq(uu)-ldosq(uu)

	 Mu=dot_product(vv(n),Mi)   
	 Mv=dot_product(vv(n+i),Mi)
	 hxy=vv(2*i+n(1))*Mp(1)+vv(2*i+n(2))*Mp(2)+vv(2*i+n(3))*Mp(3)
	 
	  do ui=1,6
  	  do uj=1,6 
	  !Posici�n 1,1 y 2,2. Cajas de las matrices de masa M afectadas por el peso SUPG. 
	  Nt(uj,ui)=Nt(uj,ui)+ccx*(Mu*Mix(uj)+Mv*Miy(uj))*Mi(ui)*jaq(uu)*jac
	  enddo								  
	  enddo	
	  
	  do uj=1,6
	  do ui=1,3
	  !Cajas de las matrices de masa N afectadas por el peso grad-div. Estas no aparecen si se usan las ecuaciones de Navier-Stokes 2D.
	  !Posicion 1,3.
	  Ntx(uj,ui)=Ntx(uj,ui)+cx*(hxy*Mix(uj)+Mhx*Mi(uj))*Mp(ui)*jaq(uu)*jac	  
	  !Posicion 2,3.
	  Nty(uj,ui)=Nty(uj,ui)+cx*(hxy*Miy(uj)+Mhy*Mi(uj))*Mp(ui)*jaq(uu)*jac	 
	  enddo
	  enddo
  
	  do ui=1,6	
	  do uj=1,3
	  !Posicion 3,3. Cajas de las matrices de masa M afectadas por el peso PSPG.	 	  	 	  
	  Ntpx(uj,ui)=Ntpx(uj,ui)+ccx*Mpx(uj)*Mi(ui)*jaq(uu)*jac
	  Ntpy(uj,ui)=Ntpy(uj,ui)+ccx*Mpy(uj)*Mi(ui)*jaq(uu)*jac 
	  enddo
	  enddo	    
	 enddo	 
	 call suma(Nt/At,n,0,0,3*i,6,6,poi)
	 call suma(Nt/At,n,i,i,3*i,6,6,podosi)	 
	 !Sin las siguientes matrices s�lo se tendr�a la estabilizaci�n SUPG con grad-div.
	 call suma(Ntx/At,n,0,2*i,3*i,6,3,potresi)  
	 call suma(Nty/At,n,i,2*i,3*i,6,3,potresi)  
	 call suma(a*Ntpx/At,n,2*i,0,3*i,3,6,poi)
	 call suma(a*Ntpy/At,n,2*i,i,3*i,3,6,podosi)
 enddo
 close(1)
deallocate(poi,podosi,potresi)
end

!-----------------------------------------------------------------------------------------------------------------------------------------
!Subrutina CAJASUPGAS
!Para las ecuaciones de aguas someras. Se calculan las matrices elementales de otros t�rminos estabilizados.
!-----------------------------------------------------------------------------------------------------------------------------------------
subroutine cajasupgas(i,j,x,y,vv,nu,del)
use elemental
integer*4, dimension(:),allocatable::poi,podosi,potresi	  
integer*4 i,j
real*8 x(i),y(i),Nn(6,6),Nm(6,6),Nk(6,3),Np(6,3),Mi(6),Mix(6),Miy(6),Mixx(6),Miyy(6),Mp(3),Mpx(3),Mpy(3)
real*8 Nh(3,6),Nf(3,6),Nl(3,3),Nj(3,6),Nw(3,6),Nxx(6,6),Nxy(6,6),Nyx(6,6),Nyy(6,6),lsic
real*8 Mu,Mv,lunoq(13),ldosq(13),jaq(13),vv(3*i),cx,ccx,nu,h,Muu,Mvv,a,b,c,hxy,Mhx,Mhy,del,Mhu 

allocate(poi(3*i),podosi(3*i),potresi(3*i))

39     format(4/,A80)
40     format(6X,6(X,I5))

!Inizializaci�n previa de variables:
!-----------------------------------
lunoq=(/0.0651301029002_8, 0.8697397941956_8, 0.0651301029002_8, 0.3128654960049_8, 0.6384441885698_8, 0.0486903154253_8, 	 &
0.6384441885698_8, 0.3128654960049_8, 0.0486903154253_8, 0.2603459660790_8, 0.4793080678419_8, 0.2603459660790_8, 1.0_8/3.0_8/)		   
ldosq=(/0.0651301029002_8, 0.0651301029002_8, 0.8697397941956_8, 0.0486903154253_8, 0.3128654960049_8, 0.6384441885698_8, 	 &
0.0486903154253_8, 0.6384441885698_8, 0.3128654960049_8, 0.2603459660790_8, 0.2603459660790_8, 0.4793080678419_8, 1.0_8/3.0_8/)
jaq=(/0.0266736178004_8, 0.0266736178004_8, 0.0266736178004_8, 0.03855688044515_8, 0.03855688044515_8, 0.03855688044515_8, 	 &
0.03855688044515_8, 0.03855688044515_8, 0.03855688044515_8, 0.0878076287166_8, 0.0878076287166_8, 0.0878076287166_8, 		 &
-0.0747850222335_8/)
do u=1,3*i
poi(u)=posi(u)
podosi(u)=posdosi(u)
potresi(u)=postresi(u)
enddo
lsic=1.0_8
a=1.0_8
  
!C�lculo de los valores para cada elemento
!-----------------------------------------
open(unit=1,file='C:\malla.txt',status='old')
read(1,39)ac
  do u=1,j
     read(1,40) n(1),n(4),n(2),n(5),n(3),n(6)

	 xa=x(n(1))-x(n(3)) 
	 xb=x(n(2))-x(n(3)) 
	 ya=y(n(1))-y(n(3))
	 yb=y(n(2))-y(n(3))  
	 jac=abs(xa*yb-xb*ya)
	 	 
	 do uu=1,6
     do uuu=1,6
     Nn(uu,uuu)=0.0_8
	 Nm(uu,uuu)=0.0_8
	 Nxx(uu,uuu)=0.0_8
	 Nxy(uu,uuu)=0.0_8
	 Nyx(uu,uuu)=0.0_8
	 Nyy(uu,uuu)=0.0_8
     enddo
	 do uuu=1,3
	 Nk(uu,uuu)=0.0_8
	 Np(uu,uuu)=0.0_8	 
     enddo
     enddo
	 do uu=1,3
     do uuu=1,6
	 Nh(uu,uuu)=0.0_8
	 Nf(uu,uuu)=0.0_8
     enddo
     enddo

	 Mixx(1)=4.0_8*(yb**2.0_8)/(jac**2.0_8)
	 Mixx(2)=4.0_8*(ya**2.0_8)/(jac**2.0_8)
	 Mixx(3)=4.0_8*((yb-ya)**2.0_8)/(jac**2.0_8)
	 Mixx(4)=-8.0_8*ya*yb/(jac**2.0_8)
	 Mixx(5)=8.0_8*(ya*yb-ya**2.0_8)/(jac**2.0_8)
	 Mixx(6)=8.0_8*(ya*yb-yb**2.0_8)/(jac**2.0_8)

	 Miyy(1)=4.0_8*(xb**2.0_8)/(jac**2.0_8)
	 Miyy(2)=4.0_8*(xa**2.0_8)/(jac**2.0_8)
	 Miyy(3)=4.0_8*((xa-xb)**2.0_8)/(jac**2.0_8)
	 Miyy(4)=-8.0_8*xa*xb/(jac**2.0_8)
	 Miyy(5)=8.0_8*(xa*xb-xa**2.0_8)/(jac**2.0_8)
	 Miyy(6)=8.0_8*(xa*xb-xb**2.0_8)/(jac**2.0_8)

	 Mpx(1)=(y(n(2))-y(n(3)))/jac
	 Mpx(2)=(y(n(3))-y(n(1)))/jac
	 Mpx(3)=(y(n(1))-y(n(2)))/jac
	 
	 Mpy(1)=(x(n(3))-x(n(2)))/jac
	 Mpy(2)=(x(n(1))-x(n(3)))/jac
	 Mpy(3)=(x(n(2))-x(n(1)))/jac

	 Mhx=vv(2*i+n(1))*Mpx(1)+vv(2*i+n(2))*Mpx(2)+vv(2*i+n(3))*Mpx(3)
	 Mhy=vv(2*i+n(1))*Mpy(1)+vv(2*i+n(2))*Mpy(2)+vv(2*i+n(3))*Mpy(3)

	 Muu=-(vv(n(1))+vv(n(2))+vv(n(3)))/9.0_8+(vv(n(4))+vv(n(5))+vv(n(6)))*4.0_8/9.0_8   
	 Mvv=-(vv(n(1)+i)+vv(n(2)+i)+vv(n(3)+i))/9.0_8+(vv(n(4)+i)+vv(n(5)+i)+vv(n(6)+i))*4.0_8/9.0_8
	 Mhu=(vv(n(1)+2*i)+vv(n(2)+2*i)+vv(n(3)+2*i))/3.0_8

	 a=sqrt((x(n(1))-x(n(2)))**2.0_8+(y(n(1))-y(n(2)))**2.0_8)
	 b=sqrt((x(n(2))-x(n(3)))**2.0_8+(y(n(2))-y(n(3)))**2.0_8)
	 c=sqrt((x(n(3))-x(n(1)))**2.0_8+(y(n(3))-y(n(1)))**2.0_8)
	 h=2.0_8*jac/(a+b+c)
 
	 cx=9.81_8*lsic/(Mhu**2.0_8)
	 ccx=9.81_8*(h**2.0_8)/(cx*(Mhu**2.0_8))  
	 
	 do uu=1,13
     
	 Mi(1)=lunoq(uu)*(2.0_8*lunoq(uu)-1.0_8)
	 Mi(2)=ldosq(uu)*(2.0_8*ldosq(uu)-1.0_8)
	 Mi(3)=(1.0_8-lunoq(uu)-ldosq(uu))*(1.0_8-2.0_8*lunoq(uu)-2.0_8*ldosq(uu))
	 Mi(4)=4.0_8*lunoq(uu)*ldosq(uu)
	 Mi(5)=4.0_8*ldosq(uu)*(1.0_8-lunoq(uu)-ldosq(uu))
	 Mi(6)=4.0_8*lunoq(uu)*(1.0_8-lunoq(uu)-ldosq(uu))

	 Mix(1)=((4.0_8*lunoq(uu)-1.0_8)*yb)/jac
	 Mix(2)=((1.0_8-4.0_8*ldosq(uu))*ya)/jac
	 Mix(3)=((4.0_8*lunoq(uu)+4.0_8*ldosq(uu)-3.0_8)*(yb-ya))/jac
	 Mix(4)=(4.0_8*ldosq(uu)*yb-4.0_8*lunoq(uu)*ya)/jac
	 Mix(5)=((8.0_8*ldosq(uu)+4.0_8*lunoq(uu)-4.0_8)*ya-4.0_8*ldosq(uu)*yb)/jac
	 Mix(6)=((4.0_8-4.0_8*ldosq(uu)-8.0_8*lunoq(uu))*yb+4.0_8*lunoq(uu)*ya)/jac
	
	 Miy(1)=((1.0_8-4.0_8*lunoq(uu))*xb)/jac
	 Miy(2)=((4.0_8*ldosq(uu)-1.0_8)*xa)/jac
	 Miy(3)=((4.0_8*lunoq(uu)+4.0_8*ldosq(uu)-3.0_8)*(xa-xb))/jac
	 Miy(4)=(4.0_8*lunoq(uu)*xa-4.0_8*ldosq(uu)*xb)/jac
	 Miy(5)=(4.0_8*ldosq(uu)*xb+(4.0_8-4.0_8*lunoq(uu)-8.0_8*ldosq(uu))*xa)/jac
	 Miy(6)=((8.0_8*lunoq(uu)+4.0_8*ldosq(uu)-4.0_8)*xb-4.0_8*lunoq(uu)*xa)/jac
 
	 Mp(1)=lunoq(uu)
	 Mp(2)=ldosq(uu)
	 Mp(3)=1.0_8-lunoq(uu)-ldosq(uu)

	 Mu=dot_product(vv(n),Mi)   
	 Mv=dot_product(vv(n+i),Mi)
	 hxy=vv(2*i+n(1))*Mp(1)+vv(2*i+n(2))*Mp(2)+vv(2*i+n(3))*Mp(3)
	 
	  do ui=1,6
  	  do uj=1,6 
	  !Posici�n 1,1 y 2,2. Cajas no lineales C con peso SUPG. 
	  Nn(uj,ui)=Nn(uj,ui)+ccx*(Mu*Mix(uj)+Mv*Miy(uj))*(Mix(ui)*Mu+Miy(ui)*Mv)*jaq(uu)*jac	  
	  !Caja A con peso SUPG (sin forma d�bil).
	  Nm(uj,ui)=Nm(uj,ui)+ccx*(Mu*Mix(uj)+Mv*Miy(uj))*(Miyy(ui)+Mixx(ui))*jaq(uu)*jac	  
	  !Cajas correspondientes a la ecuaci�n de continuidad estabilizada con pesos grad-div (se usa el par�metro 'cx'). Obviamente, son 
	  !diferentes de las usadas para las ecuaciones de Navier-Stokes 2D.
	  Nxx(uj,ui)=Nxx(uj,ui)+cx*(hxy*Mix(uj)+Mhx*Mi(uj))*(hxy*Mix(ui)+Mhx*Mi(ui))*jaq(uu)*jac	 
	  Nxy(uj,ui)=Nxy(uj,ui)+cx*(hxy*Mix(uj)+Mhx*Mi(uj))*(hxy*Miy(ui)+Mhy*Mi(ui))*jaq(uu)*jac	  
	  Nyx(uj,ui)=Nyx(uj,ui)+cx*(hxy*Miy(uj)+Mhy*Mi(uj))*(hxy*Mix(ui)+Mhx*Mi(ui))*jaq(uu)*jac	 
	  Nyy(uj,ui)=Nyy(uj,ui)+cx*(hxy*Miy(uj)+Mhy*Mi(uj))*(hxy*Miy(ui)+Mhy*Mi(ui))*jaq(uu)*jac	 
	  enddo								  
	  enddo	
	  
	  do uj=1,6
	  do ui=1,3
	  !Posicion 1,3. Caja Bx con peso SUPG.
	  !Sin forma d�bil	 	  
	  Nk(uj,ui)=Nk(uj,ui)+ccx*(Mu*Mix(uj)+Mv*Miy(uj))*Mpx(ui)*jaq(uu)*jac	  
	  !Posicion 2,3. Caja By con peso SUPG.
	  !Sin forma d�bil
	  Np(uj,ui)=Np(uj,ui)+ccx*(Mu*Mix(uj)+Mv*Miy(uj))*Mpy(ui)*jaq(uu)*jac	 
	  enddo
	  enddo
	  
	  do ui=1,6	
	  do uj=1,3
	  !Posicion 3,1. Caja C con peso PSPG.	 	  
	  Nh(uj,ui)=Nh(uj,ui)+ccx*Mpx(uj)*(Mix(ui)*Mu+Miy(ui)*Mv)*jaq(uu)*jac
	  !Posicion 3,2. Caja C con peso PSPG.	 	  
	  Nf(uj,ui)=Nf(uj,ui)+ccx*Mpy(uj)*(Mix(ui)*Mu+Miy(ui)*Mv)*jaq(uu)*jac	  	  
	  enddo
	  enddo	   
	 enddo
	 
	 !Las siguientes cajas (con pesos PSPG) no se integran ya que contienen coeficientes contantes para cada elemento	 
	 do ui=1,6	
	  do uj=1,3
	  !Posicion 3,1. Caja A con peso PSPG.	 	  
	  Nj(uj,ui)=ccx*Mpx(uj)*(Miyy(ui)+Mixx(ui))*jac/2.0_8
	  !Posicion 3,2. Caja A con peso PSPG.	 	  
	  Nw(uj,ui)=ccx*Mpy(uj)*(Miyy(ui)+Mixx(ui))*jac/2.0_8	  	  
	  enddo
	 enddo
	 do ui=1,3	
	  do uj=1,3
	  !Posicion 3,3. Caja Bx+By con peso PSPG.	 	  
	  Nl(uj,ui)=ccx*(Mpx(uj)*Mpx(ui)+Mpy(uj)*Mpy(ui))*jac/2.0_8	 	   
	  enddo
	 enddo
	 call suma((Nn-nu*Nm+Nxx)/del,n,0,0,3*i,6,6,poi)
	 call suma(Nxy/del,n,0,i,3*i,6,6,podosi)
	 call suma(Nyx/del,n,i,0,3*i,6,6,poi)
	 call suma((Nn-nu*Nm+Nyy)/del,n,i,i,3*i,6,6,podosi)
	 call suma(9.81_8*Nk/del,n,0,2*i,3*i,6,3,potresi)
	 call suma(9.81_8*Np/del,n,i,2*i,3*i,6,3,potresi)	 
	 !Sin las siguientes matrices s�lo se tendr�a la estabilizaci�n SUPG con grad-div.
	 call suma((a*Nh-a*nu*Nj)/del,n,2*i,0,3*i,3,6,poi)
	 call suma((a*Nf-a*nu*Nw)/del,n,2*i,i,3*i,3,6,podosi)
	 call suma(a*9.81_8*Nl/del,n,2*i,2*i,3*i,3,3,potresi)
 enddo
 close(1)
 deallocate(poi,podosi,potresi)
end

!-----------------------------------------------------------------------------------------------------------------------------------------
!Subrutina FSUPG
!Para las ecuaciones de aguas someras. En esta subrutina se calcula el t�rmino de fricci�n (de Manning) estabilizado.
!Se calcula para el caso yn='no' indicado en la subrutina f que eval�a las ecuaciones de aguas someras tal y como son.
!Se calculan vectores elementales que ir�n al t�rmino independiente del sistema.
!-----------------------------------------------------------------------------------------------------------------------------------------
subroutine fsupg(i,j,x,y,ma,vv,Nx,Ny,Nz) 
use elemental
integer*4 i,j	 
real*8 x(i),y(i),Nx(i),Ny(i),s,maning,ma(i),vv(3*i),Mi(6),Mu,Mv,Mui,Mvi,Mhi,dist,a,b,c
real*8 lsic,ccx,aa,Nz(i),Mix(6),Miy(6),Mpx(3),Mpy(3),h,lunoq(13),ldosq(13),jaq(13) 

39     format(4/,A80)
40     format(6X,6(X,I5))
		  													 
!Inizializaci�n previa de variables:
!-----------------------------------
lunoq=(/0.0651301029002_8, 0.8697397941956_8, 0.0651301029002_8, 0.3128654960049_8, 0.6384441885698_8, 0.0486903154253_8, 	 &
0.6384441885698_8, 0.3128654960049_8, 0.0486903154253_8, 0.2603459660790_8, 0.4793080678419_8, 0.2603459660790_8, 1.0_8/3.0_8/)		   
ldosq=(/0.0651301029002_8, 0.0651301029002_8, 0.8697397941956_8, 0.0486903154253_8, 0.3128654960049_8, 0.6384441885698_8, 	 &
0.0486903154253_8, 0.6384441885698_8, 0.3128654960049_8, 0.2603459660790_8, 0.2603459660790_8, 0.4793080678419_8, 1.0_8/3.0_8/)
jaq=(/0.0266736178004_8, 0.0266736178004_8, 0.0266736178004_8, 0.03855688044515_8, 0.03855688044515_8, 0.03855688044515_8, 	 &
0.03855688044515_8, 0.03855688044515_8, 0.03855688044515_8, 0.0878076287166_8, 0.0878076287166_8, 0.0878076287166_8, 		 &
-0.0747850222335_8/)
lsic=1.0_8
aa=1.0_8    	
	  
!C�lculo de los valores para cada elemento
!-----------------------------------------
open(unit=1,file='C:\malla.txt',status='old')
read(1,39)ac
 do u=1,j
 read(1,40) n(1),n(4),n(2),n(5),n(3),n(6)
 
	 xa=x(n(1))-x(n(3)) 
	 xb=x(n(2))-x(n(3)) 
	 ya=y(n(1))-y(n(3))
	 yb=y(n(2))-y(n(3))  
	 jac=abs(xa*yb-xb*ya)	

     !C�lculo del diametro de la circunferencia inscrita
	 a=sqrt((x(n(1))-x(n(2)))**2.0_8+(y(n(1))-y(n(2)))**2.0_8)
     b=sqrt((x(n(2))-x(n(3)))**2.0_8+(y(n(2))-y(n(3)))**2.0_8)
     c=sqrt((x(n(3))-x(n(1)))**2.0_8+(y(n(3))-y(n(1)))**2.0_8)
     h=2.0_8*jac/(a+b+c)
  
	 !C�lculo de los par�metros de estabilizaci�n para elementos que cumplen la condici�n LBB. 
	 !Aqu� s�lo es necesario calcular el par�metro 'ccx'.
     ccx=(h**2.0_8)/lsic

	 !C�lculo del m�dulo de la velocidad en cada nodo esquina y de la distancia elemental para el c�lculo del radio hidr�ulico que se estima 
	 !como el di�metro del c�rculo de �rea igual a la del elemento.
	 a=sqrt(vv(n(1))**2+vv(i+n(1))**2)
	 b=sqrt(vv(n(2))**2+vv(i+n(2))**2)
	 c=sqrt(vv(n(3))**2+vv(i+n(3))**2)
	 dist=sqrt((2.0_8*jac)/(3.14159_8))
	 
	 Mui=-(vv(n(1))+vv(n(2))+vv(n(3)))/9.0_8+(vv(n(4))+vv(n(5))+vv(n(6)))*4.0_8/9.0_8   
	 Mvi=-(vv(n(1)+i)+vv(n(2)+i)+vv(n(3)+i))/9.0_8+(vv(n(4)+i)+vv(n(5)+i)+vv(n(6)+i))*4.0_8/9.0_8
	 Mhi=(vv(n(1)+2*i)+vv(n(2)+2*i)+vv(n(3)+2*i))/3.0_8
	 maning=(ma(n(1))+ma(n(2))+ma(n(3)))/3.0_8
	 
	 if (((a.eq.0.0).and.(b.eq.0.0)).or.((a.eq.0.0).and.(c.eq.0.0)).or.((b.eq.0.0).and.(c.eq.0.0))) then
	  !Elementos pegados a un contorno con velocidades nulas. Se consideran s�lo elementos con dos nodos apoyados en la pared.
	  s=((maning**2.0_8))/(((dist*Mhi)/(dist+Mhi))**(4.0_8/3.0_8))	  
	 else 
	  !Se consideran elementos con un nodo apoyado en la pared o con ning�n nodo apoyado.
	  s=((maning**2.0_8))/(Mhi**(4.0_8/3.0_8)) 
	 endif
	 
	 Mpx(1)=(y(n(2))-y(n(3)))/jac
	 Mpx(2)=(y(n(3))-y(n(1)))/jac
	 Mpx(3)=(y(n(1))-y(n(2)))/jac
	 
	 Mpy(1)=(x(n(3))-x(n(2)))/jac
	 Mpy(2)=(x(n(1))-x(n(3)))/jac
	 Mpy(3)=(x(n(2))-x(n(1)))/jac
	  
	 do uu=1,13	      
	 
	 Mi(1)=lunoq(uu)*(2.0_8*lunoq(uu)-1.0_8)
	 Mi(2)=ldosq(uu)*(2.0_8*ldosq(uu)-1.0_8)
	 Mi(3)=(1.0_8-lunoq(uu)-ldosq(uu))*(1.0_8-2.0_8*lunoq(uu)-2.0_8*ldosq(uu))
	 Mi(4)=4.0_8*lunoq(uu)*ldosq(uu)
	 Mi(5)=4.0_8*ldosq(uu)*(1.0_8-lunoq(uu)-ldosq(uu))
	 Mi(6)=4.0_8*lunoq(uu)*(1.0_8-lunoq(uu)-ldosq(uu))

	 Mix(1)=((4.0_8*lunoq(uu)-1.0_8)*yb)/jac
	 Mix(2)=((1.0_8-4.0_8*ldosq(uu))*ya)/jac
	 Mix(3)=((4.0_8*lunoq(uu)+4.0_8*ldosq(uu)-3.0_8)*(yb-ya))/jac
	 Mix(4)=(4.0_8*ldosq(uu)*yb-4.0_8*lunoq(uu)*ya)/jac
	 Mix(5)=((8.0_8*ldosq(uu)+4.0_8*lunoq(uu)-4.0_8)*ya-4.0_8*ldosq(uu)*yb)/jac
	 Mix(6)=((4.0_8-4.0_8*ldosq(uu)-8.0_8*lunoq(uu))*yb+4.0_8*lunoq(uu)*ya)/jac
	
	 Miy(1)=((1.0_8-4.0_8*lunoq(uu))*xb)/jac
	 Miy(2)=((4.0_8*ldosq(uu)-1.0_8)*xa)/jac
	 Miy(3)=((4.0_8*lunoq(uu)+4.0_8*ldosq(uu)-3.0_8)*(xa-xb))/jac
	 Miy(4)=(4.0_8*lunoq(uu)*xa-4.0_8*ldosq(uu)*xb)/jac
	 Miy(5)=(4.0_8*ldosq(uu)*xb+(4.0_8-4.0_8*lunoq(uu)-8.0_8*ldosq(uu))*xa)/jac
	 Miy(6)=((8.0_8*lunoq(uu)+4.0_8*ldosq(uu)-4.0_8)*xb-4.0_8*lunoq(uu)*xa)/jac

	 Mu=dot_product(vv(n),Mi)   
	 Mv=dot_product(vv(n+i),Mi)
	  
	  do uuu=1,6
	  Nx(n(uuu))=Nx(n(uuu))+ccx*((Mu*Mix(uuu)+Mv*Miy(uuu))*Mu*sqrt(Mui**2+Mvi**2)*s)*jaq(uu)*jac 	 
	  Ny(n(uuu))=Ny(n(uuu))+ccx*((Mu*Mix(uuu)+Mv*Miy(uuu))*Mv*sqrt(Mui**2+Mvi**2)*s)*jaq(uu)*jac 
	  enddo
	  do uuu=1,3
	  Nz(n(uuu))=Nz(n(uuu))+aa*ccx*(Mpx(uuu)*Mu*sqrt(Mui**2+Mvi**2)*s+Mpy(uuu)*Mv*sqrt(Mui**2+Mvi**2)*s)*jaq(uu)*jac
	  enddo	  	  	
	 enddo
 enddo
 close(1)
 end

!--------------------------------------------------------------------------------------------------------------------------------------
!C�lculo de otras variables preproceso o postproceso, tras resolver las ecuaciones diferenciales. Se eval�an las funciones en los nodos
!en vez de en los puntos de integraci�n.
!--------------------------------------------------------------------------------------------------------------------------------------
!Subrutina VELOCIDADESSUBTERRANEAS
!C�lculo postproceso de velocidades	subterr�neas o velocidades de Darcy
!--------------------------------------------------------------------------------------------------------------------------------------
subroutine velocidadessubterraneas (i,j,x,y,kix,kiy,ag,vv,velx,vely)
use elemental
integer*4, dimension(:),allocatable::s
integer*4 i,j
real*8, dimension(:),allocatable::kxx,kyy,kxy
real*8 x(i),y(i),kix(i),kiy(i),vv(i),velx(i),vely(i),Mpx(3),Mpy(3),ag(i),Mpxx,Mpyy		 

allocate(s(i),kxx(i),kyy(i),kxy(i))

 39     format(4/,A80)
 40     format(6X,6(X,I5))

!Inizializaci�n previa de variables:
!-----------------------------------
do u=1,i
velx(u)=0.0_8
vely(u)=0.0_8
!Con la variable 's' se har� la media de todos los valores calculados en cada nodo.
s(u)=0
enddo
do u=1,i
kxx(u)=kix(u)*(cos(ag(u))**2)+kiy(u)*(sin(ag(u))**2)
kyy(u)=kix(u)*(sin(ag(u))**2)+kiy(u)*(cos(ag(u))**2)
kxy(u)=(kix(u)-kiy(u))*sin(ag(u))*cos(ag(u))
enddo

!C�lculo de los valores para cada nodo de cada elemento
!------------------------------------------------------
open(unit=3,file='C:\mallasub.txt',status='old')
read(3,39)ac	 
 do u=1,j
     read(3,40) n(1),n(4),n(2),n(5),n(3),n(6)
  
	 xa=x(n(1))-x(n(3)) 
	 xb=x(n(2))-x(n(3)) 
	 ya=y(n(1))-y(n(3))
	 yb=y(n(2))-y(n(3))  	 
	 jac=abs(xa*yb-xb*ya)
   	
	 Mpx(1)=(y(n(2))-y(n(3)))/jac
	 Mpx(2)=(y(n(3))-y(n(1)))/jac
	 Mpx(3)=(y(n(1))-y(n(2)))/jac
	 
	 Mpy(1)=(x(n(3))-x(n(2)))/jac
	 Mpy(2)=(x(n(1))-x(n(3)))/jac
	 Mpy(3)=(x(n(2))-x(n(1)))/jac
	 
	 Mpxx=Mpx(1)*vv(n(1))+Mpx(2)*vv(n(2))+Mpx(3)*vv(n(3))
	 Mpyy=Mpy(1)*vv(n(1))+Mpy(2)*vv(n(2))+Mpy(3)*vv(n(3))

	 !En este caso el valor de la funci�n (la derivada del nivel fre�tico) es constante en el elemento dado que el nivel fre�tico es
	 !un plano que pasa por los tres nodos del elemento y su pendiente es constante. Tendr� el mismo valor en todos los puntos.	
	 !Dado que la derivada no es cont�nua y toma distintos valores en un nodo entre distintos elementos habr� que hacer la media.	 
	 do uu=1,3
	 !Se suman valores sobre los que pueda haber y se lleva la cuenta en 's' el n�mero de coeficientes considerados por nodo.
	 velx(n(uu))=velx(n(uu))+kxx(n(uu))*Mpxx+kxy(n(uu))*Mpyy 
	 vely(n(uu))=vely(n(uu))+kxy(n(uu))*Mpxx+kyy(n(uu))*Mpyy 
	 s(n(uu))=s(n(uu))+1
	 enddo	 	 
 enddo

 !C�lculo de los valores para cada nodo
 !-------------------------------------
 do u=1,i
  !Se hace la media si se ha considerado m�s de un valor en ese nodo. En otro caso el valor ya est� bien calculado.
  if (s(u).ge.1) then
  velx(u)=-velx(u)/s(u)
  vely(u)=-vely(u)/s(u)
  endif
 enddo 
 close(3)
 deallocate(s,kxx,kyy,kxy)
end
	 
!-------------------------------------------------------------------------------------------------------------------------------
!Subrutina PENDIENTES
!C�lculo de pendientes medias del terreno (en direcciones x e y) multiplicadas por una constante 'mod' que define el usuario
!fuera de la subrutina.   !c�
!-------------------------------------------------------------------------------------------------------------------------------
subroutine pendientes (i,j,x,y,z,mod,vib)
use elemental
integer*4, dimension(:),allocatable::s
integer*4 i,j
real*8 x(i),y(i),z(i),vib(2*i),Mix(6),Miy(6),Mpxx,Mpyy,mod 		 

allocate(s(i))

 39     format(4/,A80)
 40     format(6X,6(X,I5))

!Inizializaci�n previa de variables:
!-----------------------------------
!En 'lunot,ldost', que son variables globales, van los nodos del elemento (coordenadas naturales).
!Ciertos valores son nulos dado que su dimensi�n es 7 y s�lo se tienen 6 nodos.
lunot=(/1.0_8, 0.0_8, 0.0_8, 0.5_8, 0.0_8, 0.5_8, 0.0_8/)		   
ldost=(/0.0_8, 1.0_8, 0.0_8, 0.5_8, 0.5_8, 0.0_8, 0.0_8/)
do u=1,i
vib(u)=0.0_8
vib(i+u)=0.0_8
s(u)=0
enddo

!C�lculo de los valores para cada nodo de cada elemento
!------------------------------------------------------
open(unit=1,file='C:\malla.txt',status='old')
read(1,39)ac	 
 do u=1,j
     read(1,40) n(1),n(4),n(2),n(5),n(3),n(6)
  
	 xa=x(n(1))-x(n(3)) 
	 xb=x(n(2))-x(n(3)) 
	 ya=y(n(1))-y(n(3))
	 yb=y(n(2))-y(n(3))  	 
	 jac=abs(xa*yb-xb*ya)
   	
	 do uu=1,6
	 !En este caso el valor de la funci�n (la derivada de la cota del terreno) no es constante en el elemento al evaluarse la 
	 !cota del terreno con funciones cuadr�ticas. Por tanto se eval�a en cada nodo y es necesario tener 'lunot,ldost'.  
	 Mix(1)=((4.0_8*lunot(uu)-1.0_8)*yb)/jac
	 Mix(2)=((1.0_8-4.0_8*ldost(uu))*ya)/jac
	 Mix(3)=((4.0_8*lunot(uu)+4.0_8*ldost(uu)-3.0_8)*(yb-ya))/jac
	 Mix(4)=(4.0_8*ldost(uu)*yb-4.0_8*lunot(uu)*ya)/jac
	 Mix(5)=((8.0_8*ldost(uu)+4.0_8*lunot(uu)-4.0_8)*ya-4.0_8*ldost(uu)*yb)/jac
	 Mix(6)=((4.0_8-4.0_8*ldost(uu)-8.0_8*lunot(uu))*yb+4.0_8*lunot(uu)*ya)/jac
	
	 Miy(1)=((1.0_8-4.0_8*lunot(uu))*xb)/jac
	 Miy(2)=((4.0_8*ldost(uu)-1.0_8)*xa)/jac
	 Miy(3)=((4.0_8*lunot(uu)+4.0_8*ldost(uu)-3.0_8)*(xa-xb))/jac
	 Miy(4)=(4.0_8*lunot(uu)*xa-4.0_8*ldost(uu)*xb)/jac
	 Miy(5)=(4.0_8*ldost(uu)*xb+(4.0_8-4.0_8*lunot(uu)-8.0_8*ldost(uu))*xa)/jac
	 Miy(6)=((8.0_8*lunot(uu)+4.0_8*ldost(uu)-4.0_8)*xb-4.0_8*lunot(uu)*xa)/jac
	 
	 Mpxx=dot_product(Mix,z(n))
	 Mpyy=dot_product(Miy,z(n))  
	 
	 !Dado que la derivada s�lo es cont�nua en la direcci�n de los lados del elemento habr� que hacer la media.
	 vib(n(uu))=vib(n(uu))+mod*Mpxx 
	 vib(i+n(uu))=vib(i+n(uu))+mod*Mpyy 
	 s(n(uu))=s(n(uu))+1
	 enddo	 	 
 enddo

 !C�lculo de los valores para cada nodo
 !-------------------------------------
 do u=1,i
 if (s(u).ge.1) then
 vib(u)=-vib(u)/s(u)
 vib(i+u)=-vib(i+u)/s(u)
 endif
 enddo 
 close(1)
 deallocate(s)
end

!---------------------------------------------------------------------------------------------------------------------------
!Subrutina VORTICIDAD
!C�lculo postproceso de la vorticidad (en 1/s), de las tensiones en ambas direcciones del espacio y de la tensi�n 
!tangencial(en Pa=N/m2=kg/(m*s2)).
!---------------------------------------------------------------------------------------------------------------------------
subroutine vorticidad (i,j,x,y,nu,vv,ten,vor)
use elemental
integer*4, dimension(:),allocatable::s
integer*4 i,j
real*8 x(i),y(i),vv(3*i),ten(3*i),vor(i),Mix(6),Miy(6),Mpxx,Mpyy,Mpxy,Mpyx,Mpp,nu,Mp(3)		 

allocate(s(i))

 39     format(4/,A80)
 40     format(6X,6(X,I5))

!Inizializaci�n previa de variables:
!-----------------------------------
lunot=(/1.0_8, 0.0_8, 0.0_8, 0.5_8, 0.0_8, 0.5_8, 0.0_8/)		   
ldost=(/0.0_8, 1.0_8, 0.0_8, 0.5_8, 0.5_8, 0.0_8, 0.0_8/)
do u=1,i
ten(u)=0.0_8
ten(i+u)=0.0_8
ten(2*i+u)=0.0_8
vor(u)=0.0_8
s(u)=0
enddo

!C�lculo de los valores para cada nodo de cada elemento
!------------------------------------------------------
open(unit=1,file='C:\malla.txt',status='old')
read(1,39)ac	 
 do u=1,j
     read(1,40) n(1),n(4),n(2),n(5),n(3),n(6)
  
     xa=x(n(1))-x(n(3)) 
     xb=x(n(2))-x(n(3)) 
     ya=y(n(1))-y(n(3))
     yb=y(n(2))-y(n(3))  	 
     jac=abs(xa*yb-xb*ya)
   	
     do uu=1,6

     Mp(1)=lunot(uu)
     Mp(2)=ldost(uu)
     Mp(3)=1.0_8-lunot(uu)-ldost(uu)

     Mix(1)=((4.0_8*lunot(uu)-1.0_8)*yb)/jac
     Mix(2)=((1.0_8-4.0_8*ldost(uu))*ya)/jac
     Mix(3)=((4.0_8*lunot(uu)+4.0_8*ldost(uu)-3.0_8)*(yb-ya))/jac
     Mix(4)=(4.0_8*ldost(uu)*yb-4.0_8*lunot(uu)*ya)/jac
     Mix(5)=((8.0_8*ldost(uu)+4.0_8*lunot(uu)-4.0_8)*ya-4.0_8*ldost(uu)*yb)/jac
     Mix(6)=((4.0_8-4.0_8*ldost(uu)-8.0_8*lunot(uu))*yb+4.0_8*lunot(uu)*ya)/jac
	
     Miy(1)=((1.0_8-4.0_8*lunot(uu))*xb)/jac
     Miy(2)=((4.0_8*ldost(uu)-1.0_8)*xa)/jac
     Miy(3)=((4.0_8*lunot(uu)+4.0_8*ldost(uu)-3.0_8)*(xa-xb))/jac
     Miy(4)=(4.0_8*lunot(uu)*xa-4.0_8*ldost(uu)*xb)/jac
     Miy(5)=(4.0_8*ldost(uu)*xb+(4.0_8-4.0_8*lunot(uu)-8.0_8*ldost(uu))*xa)/jac
     Miy(6)=((8.0_8*lunot(uu)+4.0_8*ldost(uu)-4.0_8)*xb-4.0_8*lunot(uu)*xa)/jac
	 
	 !C�lculo de las derivadas de la velocidad en cada nodo.
     Mpxx=dot_product(Mix,vv(n))
     Mpyy=dot_product(Miy,vv(i+n))
     Mpxy=dot_product(Miy,vv(n))
     Mpyx=dot_product(Mix,vv(i+n))  
	 !C�lculo del calado en cada nodo. S�lo habr� un valor en cada nodo (no es una derivada), pero la media se har� de todos modos 
	 !dado que no se calcular� de forma independiente el t�rmino donde el calado aparece.
	 Mpp=Mp(1)*vv(2*i+n(1))+Mp(2)*vv(2*i+n(2))+Mp(3)*vv(2*i+n(3))
	 
	 !C�lculo de tensiones en direcciones x,y y tensi�n tangencial. 
	 ten(n(uu))=ten(n(uu))+1000.0_8*(-9.81_8*Mpp+2.0_8*Mpxx*nu) 
	 ten(i+n(uu))=ten(i+n(uu))+1000.0_8*(-9.81_8*Mpp+2.0_8*Mpyy*nu)
	 ten(2*i+n(uu))=ten(2*i+n(uu))+1000.0_8*nu*(Mpxy+Mpyx) 
	 !C�lculo de la vorticidad.
	 vor(n(uu))=vor(n(uu))+Mpyx-Mpxy
	 s(n(uu))=s(n(uu))+1
	 enddo	 	 
 enddo

 !C�lculo de los valores para cada nodo
 !-------------------------------------
 do u=1,i
 if (s(u).gt.1) then
 ten(u)=ten(u)/s(u)
 ten(i+u)=ten(i+u)/s(u)
 ten(2*i+u)=ten(2*i+u)/s(u)
 vor(u)=vor(u)/s(u)
 endif
 enddo 
 close(1)
 deallocate(s)
end

!*** ---------------------------------------------------------------------------------------------------------------------------       
!*** C�lculo postproceso de velocidades que generan la conservaci�n de masa entre ambos flujos									       
!*** ---------------------------------------------------------------------------------------------------------------------------       
subroutine conservacion (i,j,x,y,z,zp,zzp,vv,velx,vely,vb)																	        !*** A�adida esta l�nea
use allocatacion																											        !*** A�adida esta l�nea
integer*4, dimension(:),allocatable::po,Esb,Ese																				        !*** A�adida esta l�nea
integer*4 i,j,u,uu,uuu,ui,w,cc,ndim,Et,a,b,c,sb(4),n(6),conv 																	        !*** A�adida esta l�nea
real*8, dimension(:),allocatable::es,esp,vec,vecdin,vbdin,vic																        !*** A�adida esta l�nea
real*8 x(i),y(i),xa,ya,z(i),zp(i),zzp(i),vv(i),velx(i),vely(i),vb(3*i),Nn(3,3),Nm(3,3),ax,ay								        !*** A�adida esta l�nea
character ac*80																												        !*** A�adida esta l�nea
logical nonzero		 																										        !*** A�adida esta l�nea
																															        !*** A�adida esta l�nea
allocate(es(i),esp(i),po(2*i),vec(2*i),Esb(2*i),Ese(2*i),vic(2*i))															        !*** A�adida esta l�nea
																															        !*** A�adida esta l�nea
 39     format(4/,A80)																										        !*** A�adida esta l�nea
 40     format(6X,6(X,I5))																									        !*** A�adida esta l�nea
																															        !*** A�adida esta l�nea
!*** Inizializaci�n previa de variables:																						       
!*** -----------------------------------																						       
sb=(/0,0,3,3/)																												        !*** A�adida esta l�nea
Et=0																														        !*** A�adida esta l�nea
do u=1,2*i																													        !*** A�adida esta l�nea
Esb(u)=0																													        !*** A�adida esta l�nea
Ese(u)=0																													        !*** A�adida esta l�nea
vec(u)=0.0_8																												        !*** A�adida esta l�nea
vic(u)=0.0_8																												        !*** A�adida esta l�nea
enddo																														        !*** A�adida esta l�nea
do u=1,i																													        !*** A�adida esta l�nea
es(u)=vv(u)-zp(u)																											        !*** A�adida esta l�nea
esp(u)=vb(2*i+u)-z(u)																											    !*** A�adida esta l�nea
enddo																															    !*** A�adida esta l�nea
 																																    !*** A�adida esta l�nea
!*** Realmente es necesario resolver un sistema para garantizar la conservaci�n.
!*** En el t�rmino independiente el c�lculo de las integrales para tres puntos con velx y espesor subt y en la matriz del sistema 
!*** ir� el calado H siendo la velocidad superficial resuelta en todos los nodos del contorno.
																																	!*** A�adida esta l�nea
!*** C�lculo del t�rmino independiente y de valores de las variables para el c�lculo de las dimensiones
!*** --------------------------------------------------------------------------------------------------
!*** S�lo se considera el contorno m�vil 
open(unit=3,file='C:\mallasub.txt',status='old')																				    !*** A�adida esta l�nea
read(3,39)ac																													    !*** A�adida esta l�nea
do u=1,j																														    !*** A�adida esta l�nea
read(3,40)n(1),n(4),n(2),n(5),n(3),n(6)																							    !*** A�adida esta l�nea
 do ui=1,3																														    !*** A�adida esta l�nea
  a=ui																															    !*** A�adida esta l�nea
  b=ui+1-sb(ui)																													    !*** A�adida esta l�nea
  c=ui+2-sb(ui+1)																												    !*** A�adida esta l�nea
  if ((zp(n(c)).ne.zzp(n(c))).and.(zp(n(b)).ne.zzp(n(b)))) then																	    !*** A�adida esta l�nea
  xa=x(n(c))-x(n(b)) 																											    !*** A�adida esta l�nea
  ya=y(n(c))-y(n(b)) 	 																										    !*** A�adida esta l�nea
  Esb(n(a))=Esb(n(a))+3																											    !*** A�adida esta l�nea
  Esb(n(c))=Esb(n(c))+3																											    !*** A�adida esta l�nea
  Esb(n(b))=Esb(n(b))+3																											    !*** A�adida esta l�nea
  Ese(n(c))=1																													    !*** A�adida esta l�nea
  Ese(n(b))=1   																												    !*** A�adida esta l�nea
  call direccion(x,y,i,n(b),n(c),n(a),ax,ay)  																					    !*** A�adida esta l�nea
  vec(n(c))=vec(n(c))+(3.0_8*velx(n(c))*es(n(c))+velx(n(b))*es(n(c))+velx(n(c))*es(n(b))+velx(n(b))*es(n(b)))*ax*abs(ya)/12.0_8	    !*** A�adida esta l�nea
  vec(n(b))=vec(n(b))+(velx(n(c))*es(n(c))+velx(n(b))*es(n(c))+velx(n(c))*es(n(b))+3.0_8*velx(n(b))*es(n(b)))*ax*abs(ya)/12.0_8     !*** A�adida esta l�nea
  vec(i+n(c))=vec(i+n(c))+(3.0_8*vely(n(c))*es(n(c))+vely(n(b))*es(n(c))+vely(n(c))*es(n(b))+vely(n(b))*es(n(b)))*ay*abs(xa)/12.0_8	!*** A�adida esta l�nea
  vec(i+n(b))=vec(i+n(b))+(vely(n(c))*es(n(c))+vely(n(b))*es(n(c))+vely(n(c))*es(n(b))+3.0_8*vely(n(b))*es(n(b)))*ay*abs(xa)/12.0_8 !*** A�adida esta l�nea
  endif																																!*** A�adida esta l�nea
 enddo																																!*** A�adida esta l�nea
enddo																																!*** A�adida esta l�nea
do u=1,i																															!*** A�adida esta l�nea
Esb(i+u)=Esb(u)																														!*** A�adida esta l�nea
Ese(i+u)=Ese(u)																														!*** A�adida esta l�nea
enddo																																!*** A�adida esta l�nea
 																																	!*** A�adida esta l�nea
!*** C�lculo de las dimensiones
!*** --------------------------
do u=1,2*i																															!*** A�adida esta l�nea
Et=Et+Esb(u)																														!*** A�adida esta l�nea
enddo																																!*** A�adida esta l�nea
Et=Et+2*i+1																															!*** A�adida esta l�nea
allocate (ita(Et),sa(Et))																											!*** A�adida esta l�nea
do u=1,Et																															!*** A�adida esta l�nea
ita(u)=0																															!*** A�adida esta l�nea
sa(u)=0.0_8																															!*** A�adida esta l�nea
enddo																																!*** A�adida esta l�nea
ita(1)=2*i+2	   																													!*** A�adida esta l�nea
do u=1,2*i																															!*** A�adida esta l�nea
ita(u+1)=ita(u)+Esb(u)																												!*** A�adida esta l�nea
enddo																																!*** A�adida esta l�nea
ndim=ita(2*i+1)-1																													!*** A�adida esta l�nea
do u=1,2*i																															!*** A�adida esta l�nea
po(u)=ita(u)																														!*** A�adida esta l�nea
enddo																																!*** A�adida esta l�nea
 																																	!*** A�adida esta l�nea
!*** C�lculo de los valores de las matrices para cada elemento 
!*** ---------------------------------------------------------
!*** S�lo se considera el contorno m�vil
rewind(3)																															!*** A�adida esta l�nea
 read(3,39)ac																														!*** A�adida esta l�nea
 do u=1,j																															!*** A�adida esta l�nea
 read(3,40)n(1),n(4),n(2),n(5),n(3),n(6)																							!*** A�adida esta l�nea
 do uu=1,3																															!*** A�adida esta l�nea
 do uuu=1,3																															!*** A�adida esta l�nea
 Nn(uu,uuu)=0.0_8																													!*** A�adida esta l�nea
 Nm(uu,uuu)=0.0_8																													!*** A�adida esta l�nea
 enddo																																!*** A�adida esta l�nea
 enddo																																!*** A�adida esta l�nea
  do ui=1,3																															!*** A�adida esta l�nea
   a=ui																																!*** A�adida esta l�nea
   b=ui+1-sb(ui)																													!*** A�adida esta l�nea
   c=ui+2-sb(ui+1)																													!*** A�adida esta l�nea
   if ((zp(n(c)).ne.zzp(n(c))).and.(zp(n(b)).ne.zzp(n(b)))) then																	!*** A�adida esta l�nea
   xa=x(n(c))-x(n(b)) 																												!*** A�adida esta l�nea
   ya=y(n(c))-y(n(b)) 																												!*** A�adida esta l�nea
   call direccion(x,y,i,n(b),n(c),n(a),ax,ay)	 																					!*** A�adida esta l�nea
   Nn(c,c)=(3.0_8*esp(n(c))+esp(n(b)))*ax*abs(ya)/12.0_8																			!*** A�adida esta l�nea
   Nn(c,b)=(esp(n(c))+esp(n(b)))*ax*abs(ya)/12.0_8																					!*** A�adida esta l�nea
   Nn(b,c)=(esp(n(c))+esp(n(b)))*ax*abs(ya)/12.0_8																					!*** A�adida esta l�nea
   Nn(b,b)=(esp(n(c))+3.0_8*esp(n(b)))*ax*abs(ya)/12.0_8   																			!*** A�adida esta l�nea
   Nm(c,c)=(3.0_8*esp(n(c))+esp(n(b)))*ay*abs(xa)/12.0_8																			!*** A�adida esta l�nea
   Nm(c,b)=(esp(n(c))+esp(n(b)))*ay*abs(xa)/12.0_8																					!*** A�adida esta l�nea
   Nm(b,c)=(esp(n(c))+esp(n(b)))*ay*abs(xa)/12.0_8																					!*** A�adida esta l�nea
   Nm(b,b)=(esp(n(c))+3.0_8*esp(n(b)))*ay*abs(xa)/12.0_8 																			!*** A�adida esta l�nea
   !*** Realmente s�lo habr� un lado que genera coeficientes (si hubiese m�s el elemento pertenece a As) y se puede poner un exit
   !*** Si hubiese m�s est� programado para tenerlo en cuenta.
   call suma(Nn,n,0,0,2*i,3,3,po)																									!*** A�adida esta l�nea
   call suma(Nm,n,i,i,2*i,3,3,po)																									!*** A�adida esta l�nea
   endif																															!*** A�adida esta l�nea
 enddo																																!*** A�adida esta l�nea
enddo																																!*** A�adida esta l�nea
																																	!*** A�adida esta l�nea
call orden (2*i,w)																													!*** A�adida esta l�nea
ndim=ndim-w																															!*** A�adida esta l�nea
																																	!*** A�adida esta l�nea
!*** C�lculo de la velocidad en contornos paralelos a x o a y
!*** --------------------------------------------------------
!*** Con el sistema cuido de que la aportaci�n sobre el vector normal de cada componente sea la misma. Sin embargo si hay contornos 
!*** paralelos a x o y una de las componentes de la velocidad no aporta valor al vector normal, y por tanto quedar�a indeterminada de 
!*** cara a la continuidad.  
!*** P.e. si tengo contorno paralelo a x queda indeterminada la velocidad en direcci�n x, podr�a ser cualquiera que no aportar� caudal y 
!*** no intervendr� en la continuidad. En este caso se tendr�a una fila de ceros en la matriz y el vector (ax nulo) y doy una velocidad 
!*** promediada considerado el cambio de espesor.  
do u=1,2*i																														    !*** A�adida esta l�nea
if ((vec(u).eq.0.0_8).and.(Ese(u).ne.0))then																						!*** A�adida esta l�nea
Ese(u)=0																															!*** A�adida esta l�nea
  if (u.le.i)then																													!*** A�adida esta l�nea
  velx(u)=velx(u)*es(u)/esp(u)																										!*** A�adida esta l�nea
  else																																!*** A�adida esta l�nea
  vely(u-i)=vely(u-i)*es(u-i)/esp(u-i)  																									!*** A�adida esta l�nea
  endif																																!*** A�adida esta l�nea
endif																																!*** A�adida esta l�nea
enddo																																!*** A�adida esta l�nea
																																	!*** A�adida esta l�nea
!*** Resoluci�n del sistema
!*** ----------------------
cc=0																																!*** A�adida esta l�nea
do u=1,2*i																															!*** A�adida esta l�nea
 if (Ese(u).ne.0) then 																												!*** A�adida esta l�nea
 !*** Tambi�n se eliminan las filas y columnas relativas a los nodos donde el contorno es paralelo a x o a y.
 cc=cc+1																															!*** A�adida esta l�nea
 vec(cc)=vec(u) 																													!*** A�adida esta l�nea
 vic(u)=sqrt(2.0_8)																													!*** A�adida esta l�nea
 endif 																																!*** A�adida esta l�nea
enddo																																!*** A�adida esta l�nea
 																																	!*** A�adida esta l�nea
call reduccion(ndim,vic,2*i)																										!*** A�adida esta l�nea
nonzero=.false.																														!*** A�adida esta l�nea
allocate(vbdin(cc),vecdin(cc))  																									!*** A�adida esta l�nea
do u=1,cc																															!*** A�adida esta l�nea
vecdin(u)=vec(u)																													!*** A�adida esta l�nea
vbdin(u)=0.0_8																														!*** A�adida esta l�nea
enddo	
conv=1																															!*** A�adida esta l�nea
call gradientesbiconjugados(vecdin,cc,vbdin,nonzero,conv)																				!*** A�adida esta l�nea
																																	!*** A�adida esta l�nea
w=0																																	!*** A�adida esta l�nea
do u=1,2*i																															!*** A�adida esta l�nea
  if (vic(u).ne.0.0_8) then																											!*** A�adida esta l�nea
  w=w+1																																!*** A�adida esta l�nea
  if (u.le.i)then																													!*** A�adida esta l�nea
  velx(u)=vbdin(w)																													!*** A�adida esta l�nea
  else																																!*** A�adida esta l�nea
  vely(u-i)=vbdin(w)																													!*** A�adida esta l�nea
  endif																																!*** A�adida esta l�nea
  endif																																!*** A�adida esta l�nea
enddo																																!*** A�adida esta l�nea
 close(3)																															!*** A�adida esta l�nea
 write(6,*)'Solucion de velocidades conservativas'																					!*** A�adida esta l�nea
 deallocate(es,esp,po,vec,Esb,Ese,vecdin,vbdin,vic,ita,sa)																			!*** A�adida esta l�nea
end																																	!*** A�adida esta l�nea
																																	!*** A�adida esta l�nea
!---------------------------------------------------------------------------------------------------------------------------------------
!Subrutinas para la construcci�n de la matriz del sistema.
!---------------------------------------------------------------------------------------------------------------------------------------
!Subrutina SUMA
!Se llama a esta subrutina cada vez que se forma una matriz elemental. Con el almacenamiento mediante vectores s�lo es necesario 
!dimensionar las matrices elementales de dimensi�n seis por seis aqu� llamadas 'Nn'. Es necesario guardar los valores enteros n(6) que 
!llevan el n�mero de los 6 nodos de cada elemento para poder colocar en su posici�n los coeficientes de cada matriz elemental. Tambi�n 
!los valores enteros 'e' y 'ee' que sit�an las cajas en el sistema en su posici�n. Aqu� sa es el vector donde van los coeficientes 
!considerados e 'ita' es el vector puntero que lleva la posici�n de estos coeficientes (se utiliza formato tipo MSR donde los 
!coeficientes no est�n ordenados por columnas, donde pueden aparecer coeficientes nulos, coeficientes de la diagonal o coeficientes de 
!la misma posici�n). 'inc' es el tama�o total del sistema creado inicialmente, antes de usar la subrutina reducci�n.
!---------------------------------------------------------------------------------------------------------------------------------------
subroutine suma(Nn,n,e,ee,inc,uj,ui,po)
use allocatacion
integer*4 n(6),r,s,rr,ss,e,ee,inc,uj,ui,po(inc)  
real*8 Nn(uj,ui)

!Escritura de coeficientes uno tras otro
!---------------------------------------
!Las primeras 'inc+1' componentes de 'ita' ya han sido referenciadas en la subrutinas dimvectas o dimvectsb.
do r=1,uj
rr=n(r)+e
do s=1,ui
ss=n(s)+ee 
sa(po(rr))=sa(po(rr))+Nn(r,s)
ita(po(rr))=ss
po(rr)=po(rr)+1
enddo
enddo
end

!------------------------------------------------------------------------------------------------------------------------
!Subrutina ORDEN
!Ordenamiento de los coeficientes almacenados y generaci�n del formato MSR para los vectores 'ita,sa'.
!Habr� que tener este formato si se quiere resolver el sistema o si se quiere hacer un producto matriz-vector.
!------------------------------------------------------------------------------------------------------------------------
subroutine orden (inc,k)
use allocatacion
integer*4 k,r,s,inc,u  
real*8 kr

!Ordenamiento de de los coeficientes de cada fila 
!------------------------------------------------
!Con este proceso aparecer�n coeficientes con la misma posici�n uno al lado del otro en el vector.
do u=1,inc
do r=ita(u),ita(u+1)-1	  
do s=r,ita(u+1)-1
 if (ita(s).lt.ita(r)) then
 k=ita(r)
 kr=sa(r)
 ita(r)=ita(s)
 sa(r)=sa(s)
 ita(s)=k
 sa(s)=kr
 endif
enddo
enddo
enddo

!Ensamblamiento y eliminaci�n de los posibles ceros y del t�rmino de la diagonal
!-------------------------------------------------------------------------------
!Se trata de procesos que cambian 'ndim'.
!Primero se ensambla sumando coeficientes con la misma posici�n. Se quitan los posibles ceros, los introducidos inicialmente 
!y los que aparecen tras el ensamblaje. Tambi�n se quita el t�rmino de la diagonal (�ste se almacena en su posici�n).  
k=0
do u=1,inc 
!En primer lugar se escribe el primer coeficiente de cada fila teniendo en cuenta la posici�n donde est�.
if ((ita(u)+k).ne.ita(u+1)) then 
!Si no hay coeficientes en la fila (si 'ita(u)' original es igual a 'ita(u+1)') no se accede a ella
ita(ita(u))=ita(ita(u)+k)
sa(ita(u))=sa(ita(u)+k)
!En segundo lugar se escriben los siguientes coeficientes de la fila ensambl�ndolos cuando sea necesario.
do r=ita(u)+k,ita(u+1)-2
 if (ita(r+1).eq.ita(r)) then
  !Se debe ensamblar el t�rmino 'r+1'.
  !Se espera a evaluar el coeficiente con el siguiente paso del bucle si se ensambla. Ser�a posible ensamblar de nuevo.	 
  sa(r-k)=sa(r-k)+sa(r+1)
  k=k+1
 else
  !No se ensambla el t�rmino 'r+1'.
  !En tercer lugar se eval�an todos los coeficientes escritos para su posible eliminaci�n.
  !Se eval�a el t�rmino 'r-k' donde se habr�n ensamblado coeficientes en caso de haber sido necesario. 
  if (ita(r-k).eq.u) then
  !Se eval�a si es un coeficiente de la diagonal.
  sa(u)=sa(r-k) 
  k=k+1
  elseif (sa(r-k).eq.0.0_8) then
  !Se eval�a si es un coeficiente nulo.
  k=k+1
  endif 
 !Se escribe el coeficiente 'r+1'. �ste no se evaluar� dentro de este bucle si se trata del �ltimo coeficiente de la fila.
 sa(r+1-k)=sa(r+1)
 ita(r+1-k)=ita(r+1)
 endif
enddo
!An�lisis del �ltimo t�rmino escrito en esa fila. Procedimiento v�lido si s�lo hay un coeficiente en la fila, en cuyo caso
!no se ha entrado en el bucle anterior pero se ha escrito inicialmente el primer y �nico coeficiente.
r=ita(u+1)-1
 if (ita(r-k).eq.u) then
 !Se eval�a si es un coeficiente de la diagonal.
 sa(u)=sa(r-k) 
 k=k+1
 elseif (sa(r-k).eq.0.0_8) then
 !Se eval�a si es un coeficiente nulo.
 k=k+1
 endif
endif
!Modificaci�n de las primeras 'inc+1' componentes del vector 'ita'.
ita(u+1)=ita(u+1)-k
enddo
end

!------------------------------------------------------------------------------------------------------------------------------------------
!Imposici�n de CC reduciendo el orden del sistema.
!------------------------------------------------------------------------------------------------------------------------------------------
!Subrutina REDUCCIONDELSISTEMA
!Se reducir� el sistema para flujo superficial (a excepci�n del que aparece con el m�todo de Newton) y para flujo subterr�neo al imponer 
!las condiciones de contorno. 
!Las condiciones de contorno vienen desde el programa principal guardadas en 'vic' si no existe un contorno m�vil. �stas habr�n sido le�das 
!del fichero 'malla.txt', del fichero 'mallasub.txt'. En otro caso las condiciones vienen desde las subrutinas que gestionan el movimiento 
!de los contornos m�viles. �stas habr�n sido le�das del fichero 'malla.txt', del fichero 'mallasub.txt' o bien dadas en ejecuci�n.
!Los contadores 'i,j' son el n�mero de nodos y de elementos que vienen desde la subrutina aguassomeras o la subrutina aguassubterranea 
!seg�n se resuelva el sistema para flujo superficial o subterr�neo. Aqu� 'inc' ser� diferente para cada sistema.
!Se modificar� la matriz almacenada en 'isa,sa' una vez que se le de formato MSR (y sea ensamblada) y el t�rmino independiente almacenado 
!en 'vector' que ya est� ensamblado. Se llamar� a la subrutina reduccion para modificar la matriz manteniendo el formato MSR.  
!-----------------------------------------------------------------------------------------------------------------------------------------
subroutine reducciondelsistema(i,c,vector,vic,v,vb,inc,bcg,ndim)
use allocatacion
integer*4, dimension(:),allocatable::ck
integer*4 i,u,uu,v(i),c,ndim,inc,w,k
real*8, dimension(:),allocatable::vi
real*8 vector(inc),vic(inc),vb(inc)
character bcg*2	   

allocate(ck(3*i),vi(inc))

!Escritura de la matriz del sistema en formato MSR (Modified Sparse Row)
!-----------------------------------------------------------------------
!Generaci�n de vectores con formato MSR.
!Se tiene que 'inc=3*i' para las ecuaciones de aguas someras o de Navier-Stokes 2D e 'inc=i' para la ecuaci�n de agua subterr�nea.
call orden (inc,k)
ndim=ndim-k

!Casos en que se reduce la matriz y el t�rmino independiente:
!Primer caso de reducci�n: Los pesos de la ecuaci�n de continuidad (en ambos modelos) s�lo est�n discretizados para los nodos esquina y 
!el calado en las din�micas (modelo superficial) o el nivel fre�tico en la ecuaci�n de continuidad (modelo subterr�neo) s�lo est�n discretizados 
!en nodos esquina. Por ello es necesario eliminar las filas y columnas de las cajas D, E, B, Asx, Asy, Ns, N,... (s�lo habr� coeficientes nulos) 
!y los respectivos coeficientes del t�rmino independiente. Se reducir� la dimensi�n del sistema de '3*i' a '2*i+nodos esquina' en el modelo 
!superficial y de 'i' a 'nodos esquina' en el subterr�neo. 
!Dar CC nulas genera el sistema que se deber�a usar con la discretizaci�n menor (elimina fila y columna simplemente).

!Segundo caso de reducci�n: Se produce al imponer las condiciones de contorno (CC) guardadas en el vector 'vic', pasando el sistema de tener 
!la dimensi�n se�alada a tener '2*i+(nodos esquina)-CC' en el modelo superficial y '(nodos esquina)-CC' en el modelo subterr�nea, donde la 
!referencia de los nodos esquina est� en el vector v(inc). Si no se usa toda la malla habr� m�s condiciones que las definidas en el contorno y 
!contorno m�vil. Esto es, la parte de la malla que no se usa con las ecuaciones superficiales (con nodos secos para el modelo conjunto o 
!superficial) o la parte de la malla que no se usa con las ecuaciones subterr�neas (modelo conjunto) no se consideran y para ello se han dado 
!CC=0 en estos nodos. 
!Las CC nulas permiten generar el sistema que se deber�a usar al no consideran las partes de la malla que tienen esa condici�n.

!Imposici�n de condiciones de contorno sobre el sistema
!------------------------------------------------------
!Se modifica el vector 'vic' de condiciones de contorno a�adiendo las CC nulas para el primer caso y as� se efect�a la reducci�n conjuntamente.
c=0
w=0
do u=1,i
 if (v(u).eq.1) then
 vic(u+inc-i)=0.0_8
 endif
enddo 

!Se modifica el t�rmino independiente con las CC almacenadas en 'vic' y se reduce el orden del t�rmino independiente (se eliminan filas).
do u=1,inc 
 if (vic(u).eq.sqrt(2.0_8)) then
 vi(u)=0.0_8
 c=c+1
 else
 vi(u)=vic(u)
 endif
enddo

do u=1,c
do while (vic(u+w).ne.sqrt(2.0_8)) 
 w=w+1
enddo
 vector(u)=vector(u+w)
 vb(u)=vb(u+w)
 do uu=ita(u+w),ita(u+w+1)-1
 vector(u)=vector(u)-sa(uu)*vi(ita(uu))
 enddo 
enddo

!Se pasa de una matriz cuadrada de dimensi�n inc a otra cuadrada con la dimensi�n se�alada en cada caso, trabajando sobre los vectores 
!'ita,sa' en que est� escrita la matriz del sistema.
call reduccion(ndim,vic,inc)
	  
!El sistema a resolver ahora es de dimension 'c' y la dimensi�n de 'ita,sa' es 'ndim' que se puede calcular como 'ita(c+1)-1' (formato MSR)
write(6,*)' '
write(6,*)'Longitud del vector matriz dispersa:',ndim

!Escritura de la matriz del sistema en formato CSC (Compressed Sparse Column)
!----------------------------------------------------------------------------
!La dimensi�n de los vectores que almacenan la matriz ser� 'ndim-1' para dos de los vectores que se puede calcular como 'ica(c+1)-1'  
!(formato CSC) y ser� 'c+1' para el tercero.			  
if (bcg.eq.'no') then
 !Se deja el formato MSR si se va a calcular con precondicionador diagonal. Los vectores 'cia,ca' no son necesarios.
 deallocate(cia,ca)
else
 !Se pasa a formato CSC si se va a calcular con precondicionador LU.
 allocate(cja(c+1))
 !Inicializaci�n de variables.
 do u=1,3*i
 ck(u)=0
 enddo
 cja(1)=1
 do u=1,c
 cja(u+1)=0
 enddo
 !C�lculo del n�mero de coeficientes por columna. 
 do u=c+2,ndim
 cja(ita(u)+1)=cja(ita(u)+1)+1
 enddo
 !Configuraci�n final de 'cja' (considerando el coeficiente de la diagonal) y escritura de la diagonal en 'cia,ca'.
 !Adem�s se tiene el contador 'ck'.
 do u=1,c						
 cja(u+1)=cja(u+1)+cja(u)+1
 ca(cja(u))=sa(u)
 cia(cja(u))=u
 ck(u)=cja(u)+1
 enddo
 !Escritura del resto de coeficientes en 'cia,ca'.
 do u=1,c										   
  do w=ita(u),ita(u+1)-1
  !En 'w' est� la columna y en 'ck(w)' la posici�n preparada para ese coeficiente.
  cia(ck(ita(w)))=u
  ca(ck(ita(w)))=sa(w)
  ck(ita(w))=ck(ita(w))+1
  enddo
 enddo
 deallocate(ita,sa)
endif
!Ahora la dimensi�n de 'ca,cia' es 'ndim-1' que se puede calcular como 'cja(c+1)-1'.
!Los vectores 'ita,sa' no son necesarios.
deallocate(ck,vi) 
end

!-----------------------------------------------------------------------------------------------------------------------------------------------
!Subrutina REDUCCI�N
!En esta subrutina se eliminan todas las filas y columnas referenciadas en 'vic' a la vez. Si existe valor en 'vic(u)' se eliminar�
!la fila y la columna 'u' del sistema existente en el momento de llamar a esta subrutina. 
!Se operar� directamente sobre los vectores 'sa,ita' que son respectivamente el vector donde van los coeficientes no nulos y el vector puntero 
!(formato MSR). 'ndim' es la dimensi�n de los vectores 'ita,sa' y 'iya,ya' guardan las primeras 'inc+1' componentes originales de 'ita,sa'.
!'lo' es el n�mero de filas-columnas a eliminar.
!'vac' es un vector que lleva en cada posici�n 'u' el n�mero de filas-columnas a eliminar con valor menor o igual que 'u' (considera 
!tambi�n la eliminaci�n de la fila-columna 'u' en ese n�mero). 
!'wo' llevar� el n�mero de coeficientes no nulos existentes en cada una de las filas a eliminar, sin considerar el coeficiente de la diagonal. 
!Incluye coeficientes correspondientes a otras columnas a eliminar.	Si la fila 'u' se debe eliminar, se deber�n eliminar 'wo(u)' coeficientes. 
!----------------------------------------------------------------------------------------------------------------------------------------------
subroutine reduccion(ndim,vic,inc)    
use allocatacion
integer*4, dimension(:),allocatable::vac,wo,iya
integer*4 u,v,w,ndim,inc,lo,le,li 
real*8, dimension(:),allocatable::ya 
real*8 vic(inc)

allocate(vac(inc),wo(inc),iya(inc+1),ya(inc+1))

!Inizializaci�n previa de variables
!----------------------------------
w=0
lo=0
do u=1,inc
 if (vic(u).ne.sqrt(2.0_8)) then
 lo=lo+1
 wo(u)=ita(u+1)-ita(u)
 else
 wo(u)=0
 endif
vac(u)=lo
enddo
do u=1,inc+1
iya(u)=ita(u)
ya(u)=sa(u)
enddo

!Se eliminan las filas 'u' del sistema si 'vic(u)' tiene valor (si es sqrt(2.0), no tiene valor)
!-----------------------------------------------------------------------------------------------
!Se sobreescriben coeficientes durante este proceso.
le=lo
do u=1,inc-lo		
 !'u' ser�n las filas a considerar (las que no ser�n eliminadas) y 'w' las que se consideran.
 do while (vic(u+w).ne.sqrt(2.0_8))
 !En 'le' (variable en el bucle) se considera la eliminaci�n los 'wo' coeficientes que se tienen hasta la fila 'u+w' (variable) y los 
 !coeficientes a eliminar en las primeras 'inc+1' componentes de 'ita,sa' originales (constante='lo').
 le=le+wo(w+u)
 w=w+1
 enddo
 !Se colocan los coeficientes (no nulos) de las filas a considerar desde la posici�n 'inc+2-lo'. Es necesario que 'le' acumule el n�mero 
 !de coeficientes almacenados en 'wo' que ya se han eliminado ya que 'le' ser� la distancia con que ser�n traspasados los coeficientes. 
 !Se escriben todos los coeficientes de las filas a considerar de forma cont�nua. De momento tendr�n todos sus coeficientes, incluyendo 
 !aqu�llos de las columnas a eliminar.  
 do v=iya(u+w)-le,iya(u+w+1)-1-le
 ita(v)=ita(v+le)
 sa(v)=sa(v+le)
 enddo
!Se va a referenciar la nueva situaci�n de los coeficientes modificando las primeras 'inc+1-lo' componentes del vector 'ita'. En ellas se 
!referenciar� la posici�n del primer coefiente (no nulo) de todas las filas a considerar. Igualmente, se tendr�n en cuenta todos sus 
!coeficientes, incluyendo aqu�llos de las columnas a eliminar. As�, 'ita(u)' siempre llevar� la posici�n inicial de la fila 'u' considerada. 
ita(u)=iya(u+w)-le
sa(u)=ya(u+w)
enddo
!Para u=inc+1-lo, 'u+w' solo es mayor que 'inc' si la �ltima fila-columna del sistema no se elimina ('w=lo'). En este caso no har� falta 
!recalcular 'le' ni 'w' y no se entrar� en el bucle. En otro caso se eliminar� la �ltima o �ltimas filas-columnas y se entrar� en el bucle.
u=inc+1-lo
do while ((u+w).le.inc)	 
le=le+wo(w+u)
w=w+1
enddo
ita(u)=iya(u+w)-le
sa(u)=0.0_8  

!Se eliminan las columnas 'u' del sistema si 'vic(u)' tiene valor 
!----------------------------------------------------------------
!'lo' continene el n�mero de filas-columnas a eliminar o coeficientes eliminados en las primeras 'inc+1' posiciones de 'ita,sa' original.
!'le' contiene el n�mero de coeficientes eliminados hasta ahora ('le-lo' son los coeficientes eliminados en las posiciones desde 'inc+2' a 
!'ndim' de 'ita,sa' originales).
!Ahora 'ita,sa' referencian a un sistema rectangular incx(inc-lo) con 'inc' columnas dado que no se han eliminado los coeficientes 
!correspondientes a las columnas a eliminar en las filas que se consideran.
!A continuaci�n se eliminan todas las columnas ('li' coeficientes) para formar un sistema (inc-lo)x(inc-lo) de forma que s�lo se pase una vez 
!por cada coeficiente de los vectores 'ita,sa'). Desde el primer coeficiente que se elimina, todos los coeficientes siguientes deben ser 
!reposicionados (a una distancia). De nuevo se sobreescriben coeficientes durante este segundo proceso.
li=0
do u=1,inc-lo										   
w=ita(u)	
 do while ((w+li).ne.ita(u+1))	  
  !Se eval�a una fila.
  if (vic(ita(w+li)).ne.sqrt(2.0_8)) then
  !Se eval�a si el coeficiente pertenece a cualquiera de las columnas a eliminar y en tal caso se modifica 'li'.
  li=li+1
  else
  !En cada fila se reposicionan los coeficientes y la referencia de su posici�n con 'vac' (sobreescribiendo sobre los coeficientes previos 
  !si pertenecen a columnas a eliminar).
  sa(w)=sa(w+li)	 	   
  ita(w)=ita(w+li)-vac(ita(w+li))		 	    
  w=w+1
  endif
 enddo
!Modificaci�n de las primeras 'inc+1-lo' componentes del vector 'ita'.
ita(u+1)=ita(u+1)-li
enddo

!C�lculo de la dimensi�n final de los vectores 'ita,sa'
!------------------------------------------------------
!'li' contiene el n�mero de coeficientes eliminados durante el segundo proceso. En total se habr�n eliminado 'le+li' coeficientes de los 
!vectores 'ita,sa' originales.
ndim=ndim-le-li
deallocate(vac,wo,iya,ya)
end

!----------------------------------------------------------------------------------------------------------------------------------------
!M�todo PBCG con precondicionamiento diagonal - Obtenida del libro Fortran recipes (chap.2). Autor:  Greenbaum, Anne, (Courant Institute)
!Precondicionador seleccionado, y c�digo reprogramado con diferentes comandos y adaptado para .f95.
!----------------------------------------------------------------------------------------------------------------------------------------
!Subrutina GRADIENTESBICONJUGADOS
!En esta subrutina se aplica m�todo PBCG con precondicionamiento diagonal (tambi�n es posible aplicar el m�todo BCG sin 
!precondicionamiento). Se resuelve el sistema lineal Ax=b. Si A es definida positiva y sim�trica se aplica el m�todo CG con m�s 
!operaciones de las necesarias.
!'i' es el n�mero de nodos de la malla, 'c' es la dimension de la matriz del sistema, 'vec' es el vector de terminos independientes del 
!sistema, 'er' indica la posibilidad de que el error est� mal calculado, 'err' es el error estimado, y 'r,s,k' son contadores. 
!Las iteraciones se paran cuando el valor 'abs(Ax-b)/abs(b)' es menor que 'tol' (tambi�n es posible utilizar otro test de convergencia
!modificando 'itol'). El n�mero m�ximo de iteraciones permitido ser� '15*c' y se toma 'tol=1e-10'. Si 'tol>err' y 'iter>15*c' se 
!obtendr� la soluci�n del sistema en cualquier caso (con mayor error).
!----------------------------------------------------------------------------------------------------------------------------------------
subroutine gradientesbiconjugados(vec,c,x,nonzero,conv)	   			   
use allocatacion
integer*4 epsil,conv
parameter(epsil=1.d-14)
integer*4 c,r,s,iter,itol,er
real*8, dimension(:),allocatable::p,pp,ra,rr,z,zz,di
real*8 err,tol,x(c),vec(c),he,ak,akden,bk,bkden,bknum,bnrm,dxnrm,xnrm,znrm,snrm 
logical nonzero	 
character tiempo*8

allocate(p(c),pp(c),ra(c),rr(c),z(c),zz(c),di(c))

302   format (A7,I8,5X,A7,E11.4E2)

!Se muestra por pantalla la hora (con formato hh:mm:ss) cuando se entra aqu�
!---------------------------------------------------------------------------
call time(tiempo)
write(6,*) 'Hora:',tiempo

!Inicializaci�n previa de variables
!----------------------------------
!Se toma 'itol=1'.
!Si 'itol=1' las iteraciones se paran cuando el valor 'abs(Ax-b)/abs(b)' es menor que 'tol'.
!Si 'itol=2' las iteraciones se paran cuando el valor 'abs(inv(P)(Ax-b))/abs(inv(P)b)' es menor que 'tol'. P es el precondicionador tal
!que P= matrizsimilar(A), y ser� mejor tanto m�s similar sea.
!Si 'itol=3' la subrutina usa su propio estimador del error en 'x', y las iteraciones paran cuando su magnitud dividida por la magnitud 
!de 'x' es menor que 'tol'. 
!Si 'itol=4' se aplica lo mismo que si 'itol=3' excepto por el hecho de que la mayor componente del error en valor absoluto y la mayor 
!componente de 'x' son utilizadas en vez de la magnitud del vector.  	  
itol=1
tol=1e-10
err=1.0_8
iter=0											   
he=0.0_8
znrm=1.0_8
dxnrm=0.0_8
xnrm=0.0_8
conv=0

!Precondicionamiento
!-------------------
!Sin precondicionamiento.
!Buenos resultados en todos los casos aunque a veces hay inestabilidades del m�todo. Se necesitan muchas iteraciones. 
!do r=1,c 
!di(r)=1.0_8     
!enddo

!Con precondicionamiento (precondicionador diagonal) o sin precondicionamiento.
!Si no hay coeficientes nulos en la diagonal de la matriz del sistema se usa la matriz diagonal P con la diagonal de dicha matriz (se 
!precondiciona con la inversa). Ocurre en GW, SW no estacionario con Picard, NS, SW si se usa Newton o estabilizaci�n.
!En el GW se mantiene o se reduce el n� iteraciones. En el SW no estacionario parece que no va bien. Para el m�todo de Newton no va bien 
!(tal vez porque los coeficientes de la diagonal en la 3� ec. son muy peque�os). Para estabilizaci�n siempre va mejor en caso estacionario 
!que si no se usa precondicionamiento. 
!En otro caso se utilizar� la matriz identidad y se resuelve sin precondicionamiento (caso anterior). Ocurre en NS, SW estacionario 
!(si no se usa Newton).
!do r=1,c
!di(r)=sa(r)
!enddo
!if (minval(abs(di)).eq.0.0_8) then    
! do r=1,c 
! di(r)=1.0_8     
! enddo 
!endif

!Con precondicionamiento (precondicionador diagonal). 
!Se usa la matriz diagonal P donde los coeficientes de la diagonal son la norma de los coeficientes de cada fila (se precondiciona con la 
!inversa). Es v�lido en todos los casos (haya o no ceros en la diagonal) como ocurre con el m�todo PCGBLU.
!Buenos resultados en todos los casos aunque a veces hay inestabilidades del m�todo. Consigue reducir el n�mero de iteraciones.
do r=1,c
di(r)=0.0_8
 do s=ita(r),ita(r+1)-1
 di(r)=di(r)+sa(s)**2
 enddo
 di(r)=sqrt(di(r)+sa(r)**2)
enddo

!Se calcula el residual inicial dependiendo de la aproximaci�n inicial dada
!--------------------------------------------------------------------------
!La aproximaci�n inicial ser� o bien un vector con sus coeficientes nulos o un vector con la �ltima soluci�n del sistema obtenida. 
  if (nonzero) then
   !Si el valor inicial 'x' no es nulo. 
   call dsprsax(x,ra,c)
   !Se realiza un producto matriz-vector con las subrutinas dsprsax.
   do s=1,c
   ra(s)=vec(s)-ra(s)
   rr(s)=ra(s)
   enddo  
  else
   !Si el valor inicial 'x' es nulo.
   do s=1,c
   ra(s)=vec(s)
   rr(s)=vec(s)
   enddo
  endif
  !Utilizando aqu� 'call dsprsax(ra,rr,c)' se calcular�a un residual 'rr' diferente y se tendr�a el algoritmo de m�nimo residual 
  !para A simetrica y no definida positiva (el GMRES para matrices no sim�tricas es de este tipo). Es una versi�n del algoritmo CG.

 !Inicializaci�n de variables seg�n el test de convergencia aplicado
 !------------------------------------------------------------------
  if (itol.eq.1) then 
   bnrm=snrm(c,vec,itol)
   !Con asolve se utiliza el precondicionador (P).
   call asolve(c,ra,z,di)
  elseif (itol.eq.2) then
   call asolve(c,vec,z,di)
   bnrm=snrm(c,z,itol)
   call asolve(c,ra,z,di)
  elseif ((itol.eq.3).or.(itol.eq.4)) then
   call asolve(c,vec,z,di)
   bnrm=snrm(c,z,itol)
   call asolve(c,ra,z,di)
   znrm=snrm(c,z,itol)
  else
   pause 'itol debe estar entre 1 y 4'
  endif 
   
   !Comienzo de las iteraciones del m�todo
   !--------------------------------------
   do while (err.gt.tol)
    !Escritura por pantalla de iteracion y del error. 
    if((iter.eq.100).or.(iter.eq.500).or.(iter.eq.c)) then   
	write(6,302) ' iter =',iter,'error =',err    
     if (er.eq.1) then
	 write(6,*)'Error posiblemente inexacto'
	 endif
	 !Parada del bucle para iter=itmax=c.
     if(iter.eq.c) then
     write(6,*)'Iteracion maxima permitida y convergencia no alcanzada'
     write(6,*)' '
	 conv=1
	 exit    
     endif 
    endif
	er=0
	iter=iter+1
    !Se debe utilizar la matriz traspuesta de la matriz precondicionamiento (PT) que es igual a la matriz de 
	!precondicionamiento (P) al ser una matriz diagonal. Por tanto se usa la subrutina asolve.
    call asolve(c,rr,zz,di)
    bknum=0.0_8
	 bknum=dot_product(z,rr)
      if (iter.eq.1) then 
       do s=1,c
       p(s)=z(s)
       pp(s)=zz(s)
       enddo
      else
       bk=bknum/bkden
       do s=1,c
       p(s)=bk*p(s)+z(s)
       pp(s)=bk*pp(s)+zz(s)
       enddo
      endif
      !A continuaci�n se calculan los residuales 'ra' y 'rr' para cada iteraci�n.
      bkden=bknum
      call dsprsax(p,z,c)
      akden=0.0_8
       akden=dot_product(z,pp)
      ak=bknum/akden
      !Se realiza un producto matriz-vector con la subrutina dsprstx.
	  call dsprstx(pp,zz,c)
	   do s=1,c
       x(s)=x(s)+ak*p(s)
       ra(s)=ra(s)-ak*z(s)
       rr(s)=rr(s)-ak*zz(s)
       enddo
	  !Con asolve se resuelve P*z=ra y se eval�a el criterio de parada.
      call asolve(c,ra,z,di)
       if (itol.eq.1) then     
        err=snrm(c,ra,itol)/bnrm
	   elseif (itol.eq.2) then
		err=snrm(c,z,itol)/bnrm
       elseif ((itol.eq.3).or.(itol.eq.4)) then
	    he=znrm
        znrm=snrm(c,z,itol)
        if (abs(he-znrm).gt.epsil*znrm) then
         dxnrm=abs(ak)*snrm(c,p,itol)
         err=znrm/abs(he-znrm)*dxnrm 
        else
	     !El error calculado podr�a ser inexacto. 
         err=znrm/bnrm
		 er=1
         cycle
        endif
         xnrm=snrm(c,x,itol)
        if (err.le.0.5_8*xnrm) then
         err=err/xnrm
        else 
	     !El error calculado podr�a ser inexacto.                        
         err=znrm/bnrm
		 er=1      
         cycle
        endif
       endif 
   enddo
   
   if (conv.ne.1) then
   write(6,*) 'Solucion obtenida en iteracion:',iter
   endif
   deallocate(p,pp,ra,rr,z,zz,di)            
end

!---------------------------------------------------------------------------------
!Funci�n SNRM
!C�lculo de la norma para el vector 'sx'.  
!---------------------------------------------------------------------------------
function snrm(c,sx,itol)
integer*4 c,itol,r,isamax
real*8 sx(c),snrm

  if (itol.le.3) then
    !C�lculo de la norma del vector.   
    snrm=sqrt(dot_product(sx,sx))
  else
    !C�lculo de la norma de la mayor componente (es directamente esa componente).
	isamax=1
    do r=1,c
     if (abs(sx(r)).gt.abs(sx(isamax))) then
     isamax=r
     endif
    enddo
    snrm=abs(sx(isamax))
  endif
end          
  
!---------------------------------------------------------------------------------------------------------------------------
!Subrutina ASOLVE
!Precondicionamiento de la matriz asociada al sistema con una matriz diagonal. Se resuelve un sistema P*x=vec con soluci�n 
!x=P(-1)*vec, estando P definido en 'di'. 
!'di' puede contener la diag(I), la diag(A) o bien la norma de los t�rminos de cada fila.
!---------------------------------------------------------------------------------------------------------------------------
subroutine asolve (c,vec,x,di)
integer*4 c,r
real*8 x(c),vec(c),di(c)   

  !Multiplico la inversa de la matriz de precondicionamiento por vec (como en la subrutina dsprsax se multiplica una matriz 
  !por x(r)). La inversa de una matriz diagonal es otra matriz diagonal con diagonal 1/(coeficientes de la diagonal). 
  !x=P(-1)*vec=vec/diag(P)=vec/di
  do r=1,c 
  x(r)=vec(r)/di(r)     
  enddo     
end
	  
!-------------------------------------------------------------------------------------------
!Subrutina DSPRSAX
!Esta subrutina multiplica una matriz dispersa de dimension 'c' almacenada por filas en 'sa' 
!con puntero 'ita' por un vector columna 'x' y lo almacena en 'vec'.
!-------------------------------------------------------------------------------------------
subroutine dsprsax(x,vec,c) 
use allocatacion
integer*4 r,k,c
real*8 vec(c),x(c)
  
  !Producto matriz-vector.
  do r=1,c
   vec(r)=sa(r)*x(r)
   do k=ita(r),ita(r+1)-1
    vec(r)=vec(r)+sa(k)*x(ita(k))
   enddo
  enddo
end  

!--------------------------------------------------------------------------------------------
!Subrutina DSPRSTX
!Esta subrutina multiplica la traspuesta de una matriz dispersa de dimension 'c' almacenada 
!por filas 'sa' con puntero 'ita' por un vector columna 'x' y lo almacena en 'vec'.
!--------------------------------------------------------------------------------------------
subroutine dsprstx(x,vec,c)
use allocatacion
integer*4 r,k,c
real*8 vec(c),x(c)
      
  !Producto matriz-vector.
  do r=1,c
   vec(r)=sa(r)*x(r) 
  enddo 
  do r=1,c 
   do k=ita(r),ita(r+1)-1     
    vec(ita(k))=vec(ita(k))+sa(k)*x(r)
   enddo
  enddo
end  

!---------------------------------------------------------------------------------------------------------------------------------------------
!M�todo PBCGLU con precondicionamiento LU - Obtenida de la librer�a SLATEC (SLAP). Autores:  Greenbaum, Anne (Courant Institute), 
!Seager, Mark K. (LLNL). C�digo reprogramado con diferentes comandos y adaptado para .f95. 
!Se definen varias constantes que dependen de la m�quina a trav�s de funciones intr�nsecas de Fortran (obviamente 1.d0=1.0_8).
!Magnitud positiva m�s peque�a = tiny(1.0). En ieee standard (intel 8087 and intel 8086 emulator) su valor es 2.22d-308.
!Magnitud m�s grande = huge(1.0). En ieee standard su valor es 1.79d308.	  
!Espacio relativo m�s peque�o = epsilon(1.0)/radix(1.0_8)=epsilon(1.0)/2.0. En ieee standard su valor es 1.11d-16.
!epsilon (1.0) ser�a el espacio relativo m�s grande.
!---------------------------------------------------------------------------------------------------------------------------------------------
!Subrutina DSLUBC
!Subrutina para resolver un sistema lineal Ax=b utilizando el m�todo de los gradientes biconjugados con descomposici�n LU incompleta 
!para formar la matriz de precondicionamiento. Los coeficientes de las matrices que se generar�n se guardar�n en 'iw y 'w' que se generan 
!aqu� y que tendr�n dimensiones 'leniw' y 'lenw' respectivamente. 'n' ser� la dimension de la matriz del sistema. 
!Las iteraciones se paran cuando el valor 'abs(Ax-b)/abs(b)' es menor que 'tol' (tambi�n es posible utilizar otro test de convergencia
!modificando 'itl'). 'itl' y 'tl' son las mismas variables que 'itol,tol' en la subrutina gradientesbiconjugados del m�todo PBCG y toman los 
!mismos valores. La subrutina devuelve 'ite' y 'er' como la iteraci�n y error estimado alcanzados. 'itmax' es el n� m�ximo de iteraciones 
!permitidas y ser� '15*c'(en la subrutina gradientesbiconjugados lo evaluaba dentro del bucle). Si 'ite>itmax' se obtendr� la soluci�n del 
!sistema en cualquier caso (con mayor error). Se obtendr� 'ite=itmax+1' si 'er>tl'.
!La subrutina devuelve 'ie' con el que se podr� saber el tipo de error que se gener� en caso de que haya alg�n problema.
!'ie=0' si todo fue bien, 'ie=1' si hay insuficiente espacio reservado para 'w,iw', 'ie=2' si no se obtiene convergencia para el n�mero m�ximo 
!de iteraciones, 'ie=3' si hay error en los datos de entrada de 'n' o 'itol', 'ie=4' si la tolerancia es demasiado estricta en cuyo caso fue
!reseteada a 500*epsilon(1.0)/2.0, 'ie=5' si la matriz de precondicionamiento no es definida positiva. Producto escalar (denominador de bk) 
!(r,z)<0. Se eval�a con una tolerancia, 'ie=6' si la matriz A no es definida positiva. Producto escalar (denominador de ak)  (p,A*p)<0. 
!Se eval�a con una tolerancia.
!Desde esta subrutina se llamar� a las subrutinas dsilus y a dbcg (en ambos casos se env�an cachos de estos vectores iw,w).
!---------------------------------------------------------------------------------------------------------------------------------------------
subroutine dslubc (b,n,x,lenw,leniw,conv)	 
use allocatacion
integer*4 crb,cib
parameter (crb=1,cib=11) 
integer*4, dimension(:),allocatable::iw
integer*4 u,ie,ite,itmax,itl,leniw,lenw,n,nt,icol,j,jbgn,jend,conv
integer*4 cdin,cdz,cil,ciu,ciw,cjl,cju,cl,cnc,cnr,cp,cpp,cr,crr,cu,cw,cz,czz,nl,nu
real*8, dimension(:),allocatable::w
real*8 er,tl,b(n),x(n) 
character tiempo*8

allocate(iw(leniw),w(lenw))
	  
!Se muestra por pantalla la hora (con formato hh:mm:ss) cuando se entra aqu�
!---------------------------------------------------------------------------
call time(tiempo)
write(6,*) 'Hora:',tiempo

!Inicializaci�n previa de variables:
!-----------------------------------
do u=1,lenw									
w(u)=0.0_8
enddo
do u=1,leniw
iw(u)=0
enddo
!Se toma 'itl=1'. Las posibilidades 'itl=1' y 'itl=2' conllevan a los mismos criterios similares a 'itol=1' y 'itol=2' en la subrutina 
!gradientesbiconjugados. 
!Si 'itl=1', las iteraciones paran cuando la norma 2 (el m�dulo) del residual dividida por la norma 2 del lado derecho es menor que 
!tol, [Ax-b]/[b]<tol. 
!Si 'itl=2', las iteraciones paran cuando la norma 2 de la inversa de M por el residual divido por la norma 2 de la inversa de M por el 
!t�rmino del lado derecho del sistema es menor que tol, [inv(M)*(Ax-b)]/[inv(M)*b]<tol.   
!M es el precondicionador tal que M=matrizsimilar(A). Sin embargo aqu� se utilizar� la matriz con la diagonal de A. 	  
itl=1
tl=1e-10
itmax=n
!nt=ndim-1 es el n�mero de elementos no nulos en A (diagonal y elem no nulos fuera de ella).
nt=cja(n+1)-1
ite=0
!Se podr�a usar la magnitud positiva m�s peque�a para inicializar er: er=tiny(1.0_8).
!Es irrelevante el valor que se defina aqu� dado que se estimar� el error antes de la primera iteraci�n en la funci�n isdbcg y se usar� 
!ese error. En gradientesbiconjugados no se estimar� antes de la primera iteraci�n y se considera que este error es 1.
er=1.0_8
ie=0
conv=0

!Con precondicionamiento (precondicionador LU). 
!Es v�lido en todos los casos (haya o no ceros en la diagonal).
!Buenos resultados en todos los casos pero parece que en mallas no regulares (backward) SW, NS con pesos BG en general hay inestabilidades 
!del m�todo para altos n�meros de Reynolds donde existir�a convergencia del sist no lineal. En este caso (Re alto) funcionar�a peor con 
!estabilizaci�n y a�n peor con Newton. Parece que en mallas regulares con SW, NS no hab�a problema con pesos BG pudiendo incluso resolverse 
!para Re alto cuando el precondicionador diagonal falla (se observ� para gran refinamiento) existiendo convergencia del sist. no lineal.   
!Con estabilizaci�n podr�a haber problemas a mayor refinamiento, y con newton con poco refinamiento (Re alto).
!Se tiene un menor n�mero de iteraciones que con m�todo PCGB con precondicionador diagonal.
	 	  
!Comprueba si hay sistema (siempre habr�)
!----------------------------------------
if ((n.lt.1).or.(nt.lt.1)) then
ie=3
write(6,*)'Se tiene el error (buscar por numero):',ie
return
endif

 !Cuenta el n�mero de coeficientes no nulos de la matriz de precondicionamiento ilu.
 !---------------------------------------------------------------------------------- 
 !C�lculo de 'nl,nu'.
  nl=0
  nu=0
  do icol=1,n						 
    !No cuenta en la diagonal.
    jbgn=cja(icol)+1
    jend=cja(icol+1)-1
    if(jbgn.le.jend) then
    !Si hay coeficientes no nulos (diferentes del de la diagonal) en la columna (A est� almacenada por columnas).
      do j=jbgn,jend
        if(cia(j).gt.icol) then
		 !Coeficientes de la matriz triangular inferior.
         nl=nl+1
        else
		 !Coeficientes de la matriz triangular superior.
         nu=nu+1
        endif
      enddo
    endif
  enddo

 !A continuaci�n establecen las matrices de trabajo.
 !--------------------------------------------------
 cil=cib
 cjl=cil+n+1
 ciu=cjl+nl
 cju=ciu+nu
 cnr=cju+n+1
 cnc=cnr+n
 ciw=cnc+n

 cl=crb
 cdin=cl+nl
 cu=cdin+n
 cr=cu+nu
 cz=cr+n
 cp=cz+n
 crr=cp+n
 czz=crr+n
 cpp=czz+n
 cdz=cpp+n
 cw=cdz+n

 !Chequea el espacio reservado para los vectores 'w,iw'.
 call dchkw(ciw,leniw,cw,lenw,ie,ite) 
 if (ie.eq.1) then
 write(6,*)'Se tiene el error (buscar por numero):',ie
 return
 endif

 !Si a�n no hay error previos se referencia en el vector puntero 'iw' en qu� partes de los vectores 'iw,w' ir�n
 !las matrices que se almacenar�n.
 iw(1)=cil
 iw(2)=cjl
 iw(3)=ciu
 iw(4)=cju
 iw(5)=cl
 iw(6)=cdin
 iw(7)=cu
 iw(9)=ciw
 iw(10)=cw

 !A continuaci�n se computa la descomposici�n lu incompleta.
 !----------------------------------------------------------
 !Sea A, de 10 componentes. A(2:4) es un vector con la primera componente en 2 y la �ltima en 4 de modo que tiene 3 componentes.
 !Cuando se pasa a una subrutina A(2) se env�a un trozo del vector A. El vector que se genera en la siguiente tendr� una dimensi�n 8.
 !As� pues, se env�an a la subrutina dsilus trozos de los vectores.
 call dsilus(n,nl,iw(cil),iw(cjl),w(cl),w(cdin),nu,iw(ciu),iw(cju),w(cu),iw(cnr),iw(cnc))	           																                                               
 
 !A continuaci�n se efect�a el algoritmo de los gradientes biconjugados con el precondicionamiento indicado.
 !----------------------------------------------------------------------------------------------------------
 !Se env�a a la subrutina dbcg otros trozos de los vectores. Desde ella se llamar� de nuevo a varias subrutinas.
 call dbcg(n,b,x,itl,tl,itmax,ite,er,ie,w(cr),w(cz),w(cp),w(crr),w(czz),w(cpp),w(cdz),w,iw,lenw,leniw,conv)
  
 !Se muestra por pantalla el error obtenido en caso de que haya alg�n problema
 !---------------------------------------------------------------------------- 
 if (ie.eq.0) then
 write(6,*) 'Solucion obtenida en iteracion:',ite
 else
 if (ie.eq.2) then
 write(6,*)'Iteracion maxima permitida y convergencia no alcanzada'
 write(6,*)' '
 else
 write(6,*)'Se tiene el error (buscar por numero):',ie
  if (ite.gt.0) then
  write(6,*)'error=',er
  endif
 endif
 endif
deallocate(iw,w)
end

!----------------------------------------------------------------------------------------------------------------------------------
!Subrutina DCHKW
!Esta subrutina chequea si la longitud reservada para las matrices de trabajo (iw,w) es suficiente y devuelve el error 
!correspondiente en otro caso. 
!----------------------------------------------------------------------------------------------------------------------------------
subroutine dchkw (lociw,leniw,locw,lenw,ierr,iter) 
integer*4 ierr,iter,leniw,lenw,lociw,locw

!Inicializaci�n previa de variables
!----------------------------------
ierr=0
iter=0
  
!Chequeo de la longitud reservada
!--------------------------------
!Chequea el espacio de la matriz de coeficientes enteros de trabajo.
if(lociw.gt.leniw) then 	 
   ierr=1
   write(6,*) 'Espacio necesario matriz entera:',lociw
   write(6,*) 'Espacio disponible:',leniw
endif
!Chequea el espacio de la matriz de coeficientes reales de trabajo.
if(locw.gt.lenw) then 
   ierr=1
   write(6,*) 'Espacio necesario matriz real:',locw
   write(6,*) 'Espacio disponible:',lenw
endif
end

!--------------------------------------------------------------------------------------------------------------------
!Subrutina DSLUI
!Esta subrutina act�a como una interfaz entre una subrutina y la subrutina que realmente computa x=inv(LDU)b.
!iw(1)=localizaci�n del primer coeficiente de 'il' en iw, iw(2)=localizaci�n del primer coeficiente de 'jl' en iw.
!iw(3)=localizaci�n del primer coeficiente de 'iu' en iw, iw(4)=localizaci�n del primer coeficiente de 'ju' en iw.
!iw(5)=localizaci�n del primer coeficiente de 'l' en rw,  iw(6)=localizaci�n del primer coeficiente de 'dinv' en rw.
!iw(7)=localizaci�n del primer coeficiente de 'u' en rw.
!--------------------------------------------------------------------------------------------------------------------
subroutine dslui (n,b,x,rw,iw,lw,liw)
use allocatacion
integer*4 n,lw,liw,nu,nl,iw(liw),locdin,locil,lociu,locjl,locju,locl,locu
real*8 b(n),rw(lw),x(n) 

!Saca la ubicaci�n donde se escribir�n las matrices que sostienen la factorizaci�n ILU.
!--------------------------------------------------------------------------------------
locil=iw(1)
locjl=iw(2)
lociu=iw(3)
locju=iw(4)
locl=iw(5)
locdin=iw(6)
locu=iw(7)
nl=iw(3)-iw(2)
nu=iw(4)-iw(3)

!Resuelve el sistema LUx=b
!-------------------------
call dslui2(n,b,x,nu,nl,iw(locil),iw(locjl),rw(locl),rw(locdin),iw(lociu),iw(locju),rw(locu))
end

!---------------------------------------------------------------------------------------------------------------------------
!Subrutina DSLUI2
!Subrutina para resolver un sistema de la forma L*D*Ux=b, donde L es una matriz triangular inferior, D es una matriz 
!diagonal y U es una matriz triangular superior.
!'il,jl,l' contendr�n la unidad triangular inferior L de la descomposici�n incompleta de la matriz A, almacenada 
!por filas (formato CSR similar al CSC). 'iu,ju,u' contendr�n la unidad triangular superior U de la descomposici�n 
!incompleta de la matriz A, almacenada por columnas (formato CSC). Esta factorizaci�n ILU es computada por la subrutina 
!dsilus. La diagonal (todos sus coeficientes tendr�n valor 1) es almacenada.
!--------------------------------------------------------------------------------------------------------------------------- 
subroutine dslui2 (n,b,x,nu,nl,il,jl,l,dinv,iu,ju,u)
integer*4 n,nu,nl,il(n+1),iu(nu),jl(nl),ju(n+1),i,icol,irow,j,jbgn,jend
real*8 b(n),dinv(n),x(n),l(nl),u(nu) 

!Se almacena el t�rmino independiente (b) en 'x'
!-----------------------------------------------
do i=1,n
  x(i)=b(i)
enddo

!Resuelve L*y=b, almacenando el resultado (y) en 'x' (L almacenada por filas).
!-----------------------------------------------------------------------------
do irow=2,n
  jbgn=il(irow)
  jend=il(irow+1)-1
  if(jbgn.le.jend) then
    do j=jbgn,jend
      x(irow)=x(irow)-l(j)*x(jl(j))
    enddo
  endif
enddo

!Resuelve D*z=y, almacenando el resultado (z) en 'x'.
!----------------------------------------------------
do i=1,n
  x(i)=x(i)*dinv(i)
enddo

!Resuelve U*x=z, (U almacenada por columnas).
!--------------------------------------------	  
do icol=n,2,-1
  jbgn=ju(icol)
  jend=ju(icol+1)-1
  if(jbgn.le.jend) then
    do j=jbgn,jend
      x(iu(j))=x(iu(j))-u(j)*x(icol)
    enddo
  endif
enddo
end

!------------------------------------------------------------------------------------------------------------------------------
!Subrutina DSILUS
!Subrutina para generar la descomposici�n incompleta LDU de la matriz siendo L la unidad triangular inferior, U la unidad
!triangular superior. Se almacena la inversa de la diangonal D.	 
!'nl' es el n�mero de coeficientes no nulos en la matriz L, 'nu' es el n�mero de coeficientes no nulos en la matriz U.
!'nrow(i)' es el n�mero de coeficientes no nulos en la fila 'i' de L, 'ncol(i)' es el n�mero de coeficientes no nulos 
!en la columna 'i' de U.
!------------------------------------------------------------------------------------------------------------------------------ 
subroutine dsilus (n,nl,il,jl,l,dinv,nu,iu,ju,u,nrow,ncol)
use allocatacion
integer*4 n,nl,nu,il(n+1),iu(nu),jl(nl),ju(n+1),ncol(n),nrow(n),i,ibgn,icol,iend,indx,indx1,indx2
integer*4 indxc1,indxc2,indxr1,indxr2,irow,itemp,j,jbgn,jend,jtemp,k,kc,kr,kk 
real*8 dinv(n),l(nl),u(nu),temp
      
!El significado de las variables (cu�l lleva la posici�n y cu�l lleva la informaci�n de donde empieza cada columna) enteras 
!que llevan la matriz "jmatriz y imatriz" cambia entre la matriz triangular inferior y superior.
!'cja,cia,ca' (formato columna) equivale a 'ju,iu,u' (formato columna) y equivale a 'il,jl,l' (formato fila). 

  !C�lculo del n�mero de coeficientes en cada fila de L y en cada columna de U.
  !----------------------------------------------------------------------------
   do i=1,n
     nrow(i)=0
     ncol(i)=0
   enddo
   do icol=1,n
     jbgn=cja(icol)+1
     jend=cja(icol+1)-1
     if(jbgn.le.jend) then
	 !Si hay coeficientes en esa columna.
       do j=jbgn,jend
         if(cia(j).lt.icol) then
          !C�lculo del n�mero de elementos no nulos de la matriz triangular superior en la columna 'icol'.
		  !Se a�aden todos para cada 'icol'.
          ncol(icol)=ncol(icol)+1
         else
	      !cia(j) contiene la fila en que est� posicionado el coeficiente dentro de la columna 'icol'.
		  !C�lculo del n�mero de elementos no nulos en la columna de la matriz triangular inferior.
	      !Se van a�adiendo elementos en la misma fila para varios 'icol'.
          nrow(cia(j))=nrow(cia(j))+1
         endif
       enddo
     endif
   enddo
   
   !Copia la matriz A dentro de las estructuras L y U.
   !--------------------------------------------------
   !Se arreglan los vectores punteros 'ju,il'
   ju(1)=1
   il(1)=1
   do icol=1,n
     il(icol+1)=il(icol)+nrow(icol)
     ju(icol+1)=ju(icol)+ncol(icol)
     nrow(icol)=il(icol)
     ncol(icol)=ju(icol)
   enddo
   !Se escriben los coeficientes y se arreglan los vectores de componentes enteras 'iu,jl'.
   do icol=1,n
     dinv(icol)=ca(cja(icol))
     jbgn=cja(icol)+1
     jend=cja(icol+1)-1
     if(jbgn.le.jend) then
       do j=jbgn,jend
         irow=cia(j)
         if(irow.lt.icol) then
           !Coeficientes del tri�ngulo superior.
           iu(ncol(icol))=irow
           u(ncol(icol))=ca(j)
           ncol(icol)=ncol(icol)+1
         else
           !Coeficientes del tri�ngulo inferior (almacenado por filas).
           jl(nrow(irow))=icol
           l(nrow(irow))=ca(j)
           nrow(irow)=nrow(irow)+1
         endif
       enddo
     endif
   enddo
   !L est� almacenada en las primeras 'nl' componentes de 'l,jl,il' y U est� almacenada en las primeras 'nu' 
   !componentes de 'u,iu,ju'.

   !Se ordenan las filas de L y las columnas de U (creo que ya est�n ordenadas y este procedimiento no es necesario)
   !----------------------------------------------------------------------------------------------------------------
   do k=2,n
     jbgn=ju(k)
     jend=ju(k+1)-1
     if(jbgn.lt.jend) then
       do j=jbgn,jend-1
         do i=j+1,jend
           if(iu(j).gt.iu(i)) then
             itemp=iu(j)
             iu(j)=iu(i)
             iu(i)=itemp
             temp=u(j)
             u(j)=u(i)
             u(i)=temp
           endif
         enddo
       enddo
     endif
     ibgn=il(k)
     iend=il(k+1)-1
     if(ibgn.lt.iend) then
       do i=ibgn,iend-1
         do j=i+1,iend
           if(jl(i).gt.jl(j)) then
             jtemp=ju(i)
             ju(i)=ju(j)
             ju(j)=jtemp
             temp=l(i)
             l(i)=l(j)
             l(j)=temp
           endif
         enddo
       enddo
     endif
   enddo

   !Se realiza la descomposici�n LDU incompleta.
   !--------------------------------------------
   do i=2,n
   !Se eval�a la fila 'i' de L.
     indx1=il(i)
     indx2=il(i+1)-1
     if(indx1.le.indx2) then
     do indx=indx1,indx2
       if(indx.ne.indx1) then
       indxr1=indx1
       indxr2=indx-1
       indxc1=ju(jl(indx))
       indxc2=ju(jl(indx)+1)-1
       if(indxc1.le.indxc2) then
 		  kk=1	
		  kr=jl(indxr1)
		  do while (kk.eq.1)	
		  kk=0
		  kc=iu(indxc1) 
		   if(kr.gt.kc) then
           indxc1=indxc1+1
            if(indxc1.le.indxc2) then     
			kk=1
			endif	      
		   elseif(kr.lt.kc) then
           indxr1=indxr1+1
            if(indxr1.le.indxr2) then
			kr=jl(indxr1)
			kk=1
			endif
           else
           l(indx)=l(indx)-l(indxr1)*dinv(kc)*u(indxc1)
           indxr1=indxr1+1 
           indxc1=indxc1+1 
            if((indxr1.le.indxr2).and.(indxc1.le.indxc2)) then
			kr=jl(indxr1)
			kk=1
			endif
		   endif
 		  enddo
	   endif
	   endif
       l(indx)=l(indx)/dinv(jl(indx))
     enddo
	 endif

     !Se eval�a la columna 'i' de U.
     indx1=ju(i)
     indx2=ju(i+1)-1
     if(indx1.le.indx2) then
     do indx=indx1,indx2
       if(indx.ne.indx1) then
       indxc1=indx1
       indxc2=indx-1	
       indxr1=il(iu(indx))
       indxr2=il(iu(indx)+1)-1
       if(indxr1.le.indxr2) then
		  kk=1
		  kr=jl(indxr1)
		  do while (kk.eq.1)	
		  kk=0
		  kc=iu(indxc1)
		   if(kr.gt.kc) then
           indxc1=indxc1+1
            if(indxc1.le.indxc2) then
			kk=1
			endif
           elseif(kr.lt.kc) then
           indxr1=indxr1+1
            if(indxr1.le.indxr2) then			   
			kr=jl(indxr1)
			kk=1
			endif
           else
           u(indx)=u(indx)-l(indxr1)*dinv(kc)*u(indxc1)
           indxr1=indxr1+1
           indxc1=indxc1+1
            if((indxr1.le.indxr2).and.(indxc1.le.indxc2)) then			   
			kr=jl(indxr1)
			kk=1
			endif
           endif
		  enddo
	   endif
	   endif
       u(indx)=u(indx)/dinv(iu(indx))
     enddo
	 endif

     !Se eval�a el coeficiente 'i' de la diagonal.
     indxr1=il(i)
     indxr2=il(i+1)-1
     if(indxr1.le.indxr2) then  
     indxc1=ju(i)
     indxc2=ju(i+1)-1
     if(indxc1.le.indxc2) then 
	   kk=1
	   kr=jl(indxr1)
	   do while (kk.eq.1)	
	   kk=0
	   kc=iu(indxc1)
		if(kr.gt.kc) then
        indxc1=indxc1+1
         if(indxc1.le.indxc2) then
		 kk=1
		 endif
        elseif(kr.lt.kc) then
        indxr1=indxr1+1
         if(indxr1.le.indxr2) then			   
		 kr=jl(indxr1)
		 kk=1
		 endif
        else
        dinv(i)=dinv(i)-l(indxr1)*dinv(kc)*u(indxc1)
        indxr1=indxr1+1
        indxc1=indxc1+1
         if((indxr1.le.indxr2).and.(indxc1.le.indxc2)) then			
		 kr=jl(indxr1)
		 kk=1
		 endif
        endif
	   enddo  
	 endif
	 endif
   enddo

   !Se reemplazan los elementos diagonales por sus inversos. 
   do i=1,n
    dinv(i)=1.0_8/dinv(i)
   enddo
end

!-------------------------------------------------------------------------------------------------------------------------
!Subrutina DSLUTI
!Esta subrutina act�a como una interfaz entre una subrutina y la subrutina que realmente computa x=inv(transpuesta(LDU))b.
!iw(1)=localizaci�n del primer coeficiente de 'il' en 'iw', iw(2)=localizaci�n del primer coeficiente de 'jl' en 'iw'.
!iw(3)=localizaci�n del primer coeficiente de 'iu' en 'iw', iw(4)=localizaci�n del primer coeficiente de 'ju' en 'iw'.
!iw(5)=localizaci�n del primer coeficiente de 'l' en 'rw',  iw(6)=localizaci�n del primer coeficiente de 'dinv' en 'rw'.
!iw(7)=localizaci�n del primer coeficiente de 'u' en 'rw'.
!-------------------------------------------------------------------------------------------------------------------------
subroutine dsluti (n,b,x,rw,iw,lw,liw)
integer*4 n,lw,liw,nu,nl,iw(liw),locdin,locil,lociu,locjl,locju,locl,locu
real*8 b(n),rw(lw),x(n) 

!Saca los punteros a las matrices L,D y U.
!-----------------------------------------
locil=iw(1)
locjl=iw(2)
lociu=iw(3)
locju=iw(4)
locl=iw(5)
locdin=iw(6)
locu=iw(7)
nl=iw(3)-iw(2)
nu=iw(4)-iw(3)

!Resuelve el sistema transpuesta(LDU)x=b
!---------------------------------------
!Los arrays con () dentro de un call son nuevos arrays cuya primera componente es la que va entre par�ntesis.  
call dslui4(n,b,x,nu,nl,iw(locil),iw(locjl),rw(locl),rw(locdin),iw(lociu),iw(locju),rw(locu))
end

!-------------------------------------------------------------------------------------------------------------------------
!Subrutina DSLUI4
!Subrutina para resolver un sistema de la forma transpuesta(L*D*U)x=b, donde L es una matriz triangular inferior 
!D es una matriz diagonal y U es una matriz triangular superior.
!-------------------------------------------------------------------------------------------------------------------
subroutine dslui4 (n,b,x,nu,nl,il,jl,l,dinv,iu,ju,u)							   
integer*4 n,nu,nl,il(n+1),iu(nu),jl(nl),ju(n+1),i,icol,irow,j,jbgn,jend
real*8 b(n),dinv(n),x(n),l(nl),u(nu) 

!Se almacena el t�rmino independiente (b) en 'x'
!-----------------------------------------------
do i=1,n
  x(i)=b(i)
enddo

!Resuelve transpuesta(U)*y=x, almacenando el resultado (y) en 'x' (U almacenada por columnas).	
!---------------------------------------------------------------------------------------------
do irow=2,n
  jbgn=ju(irow)
  jend=ju(irow+1)-1
  if(jbgn.le.jend) then
    do j=jbgn,jend
      x(irow)=x(irow)-u(j)*x(iu(j))
    enddo
  endif
enddo

!Resuelve D*z=y, almacenando el resultado (z) en 'x'.	
!---------------------------------------------------- 
!'dinv' es la inversa de la matriz diagonal D.
do i=1,n
  x(i)=x(i)*dinv(i)
enddo

!Resuelve transpuesta(L)*x=z (L almacenada por filas).
!-----------------------------------------------------
do icol=n,2,-1
  jbgn=il(icol)
  jend=il(icol+1)-1
  if(jbgn.le.jend) then
    do j=jbgn,jend
      x(jl(j))=x(jl(j))-l(j)*x(icol)
    enddo
  endif
enddo
end

!--------------------------------------------------------------------------------------------------------------------------------
!Subrutina DSMTV
!Subrutina para calcular el producto matriz-vector: y = transpuesta(A)*x. La matriz dispersa A estar� almacenada en formato CSC.
!Es equivalente a la subrutina dsprstx (considera formato MSR para A). 
!--------------------------------------------------------------------------------------------------------------------------------
subroutine dsmtv (n,x,y) 	 
use allocatacion
integer*4 n,i,irow
real*8 x(n),y(n) 

!Inicializa a cero el vector resultado.
!--------------------------------------
do i=1,n
  y(i)=0.0_8
enddo

!Multiplica por transpuesta(A).
!------------------------------ 
!Se tiene en cuenta que se tiene transpuesta(A) si se considera que A est� almacenada por filas.
do irow=1,n
  do i=cja(irow),cja(irow+1)-1
    y(irow)=y(irow)+ca(i)*x(cia(i))
  enddo
enddo
end

!---------------------------------------------------------------------------------------------------------------------------------
!Subrutina DSMV
!Subrutina para calcular el producto matriz-vector: y = A*x. Es	equivalente a la subrutina dsprsax (considera formato MSR para A).
!---------------------------------------------------------------------------------------------------------------------------------
subroutine dsmv (n,x,y)
use allocatacion
integer*4 n,i,icol
real*8 x(n),y(n) 

!Inicializa a cero el vector resultado.
!--------------------------------------
do i=1,n
  y(i)=0.0_8
enddo

!Multiplica por A.
!-----------------
do icol=1,n
  do i=cja(icol),cja(icol+1)-1
    y(cia(i))=y(cia(i))+ca(i)*x(icol)
  enddo
enddo
end

!---------------------------------------------------------------------------------------------------------------------------
!Subrutina DBCG
!Subrutina para resolver un sistema lineal no sim�trico Ax=b usando el m�todo de los gradientes biconjugados precondicionado
!Se llevan a cabo las iteraciones del m�todo. Es equivalente a la subrutina gradientesbiconjugados.
!'w' contiene por este orden: (l,dinv,u,r,z,p,rr,zz,pp,dz) con dimensiones (nl,n,nu,n,n,n,n,n,n,n) respectivamente.
!'iw' contiene por este orden: (�ndices,il,jl,iu,ju,nrow,ncol) con dimensiones (10,n+1,nl,nu,n+1,n) respectivamente.
!---------------------------------------------------------------------------------------------------------------------------
subroutine dbcg(n,b,x,itol,tol,itmax,iter,err,ierr,r,z,p,rr,zz,pp,dz,rw,iw,lenw,leniw,conv)   
use allocatacion
integer*4 i,k,ierr,iter,itmax,itol,n,lenw,leniw,iw(leniw),isdbcg,conv
real*8 err,tol,b(n),dz(n),p(n),pp(n),r(n),rr(n),rw(lenw),x(n),z(n),zz(n)  
real*8 ak,akden,bk,bkden,bknum,bnrm,tolmin,fuzz   			   

!Inicializaci�n previa de variables
!----------------------------------
iter=0
!Se usa el espacio relativo m�s peque�o  
fuzz=epsilon(1.0_8)/2.0_8	
tolmin=500.0_8*fuzz
fuzz=fuzz*fuzz
!Chequeo de datos de entrada:
if(tol.lt.tolmin) then
  tol=tolmin
  ierr=4
endif

 !Se calcula el residual inicial y el pseudo-residual 
 !---------------------------------------------------  
 call dsmv(n,x,r)
 do i=1,n
    r(i)=b(i)-r(i)
    rr(i)=r(i)
 enddo
 call dslui(n,r,z,rw,iw,lenw,leniw)
 call dsluti(n,rr,zz,rw,iw,lenw,leniw)

 !Si se cumple el criterio de parada se sale de la subrutina sin hacer iteraci�n alguna.
 if(isdbcg(n,b,itol,tol,iter,err,ierr,r,z,dz,rw,iw,bnrm,lenw,leniw).ne.0) then
 return
 endif
 if(ierr.ne.0) then
 return
 endif

 !Comienzo de las iteraciones del m�todo
 !--------------------------------------
 do k=1,itmax
   iter=k

   !Calcula el coeficiente 'bk' y los vectores de direcci�n 'p' y 'pp'.
   bknum=0.0_8
   !Se calcula el producto interior de dos vectores.
   bknum=dot_product(z,rr)
   if(abs(bknum).le.fuzz) then
     ierr=5
     return
   endif
   if(iter.eq.1) then
     !Se hace la copia un vector.
     do i=1,n
     p(i)=z(i)
     pp(i)=zz(i)
     enddo
   else
     bk=bknum/bkden
     do i=1,n
       p(i)=z(i)+bk*p(i)
       pp(i)=zz(i)+bk*pp(i)
     enddo
   endif
   bkden=bknum

   !Calcula el coeficiente 'ak', nueva iteraci�n de 'x', nuevos residuales 'r' y 'rr' y 
   !nuevos pseudo-residuales 'z' y 'zz'.
   call dsmv(n,p,z)
   akden=0.0_8
   !Se calcula el producto interior de dos vectores.
   akden=dot_product(pp,z)	 
   ak=bknum/akden
   if(abs(akden).le.fuzz) then
     ierr=6
     return
   endif
   !Se calcula la suma de una constante 'ak' multiplicada por un vector y un vector.
   if (ak.ne.0.0_8) then
    do i=1,n
    x(i)=x(i)+ak*p(i)
    r(i)=r(i)-ak*z(i)
    enddo
   endif
   call dsmtv(n,pp,zz)
   !Se calcula la suma de una constante (ak) multiplicada por un vector y un vector
   if (ak.ne.0.0_8) then
    do i=1,n
    rr(i)=rr(i)-ak*zz(i)
    enddo
   endif
   call dslui(n,r,z,rw,iw,lenw,leniw)
   call dsluti(n,rr,zz,rw,iw,lenw,leniw)

   !Eval�a el criterio de parada.
   if (isdbcg(n,b,itol,tol,iter,err,ierr,r,z,dz,rw,iw,bnrm,lenw,leniw).ne.0) then 
   return
   endif
 enddo

!Criterio de parada no satisfecho.
iter=itmax+1
ierr=2
conv=1
end   

!-------------------------------------------------------------------------------------------------------------------------------------
!Funci�n ISDBCG
!Esta subrutina eval�a la parada para el esquema iterativo BCG. Devuelve un valor no nulo si el error estimado (tipo de 
!estimaci�n determinada por 'itol') es menor que la tolerancia especificada (tol).
!'bnrm' es la norma del t�rmino del lado derecho del sistema considerado (depende de 'itol'). Es calculado s�lo en la primera llamada. 
!Valores que devuelve la funci�n:
! 0 : El estimador del error no es menor que la tolerancia especificada. Las iteraciones deben continuar.
! 1 : El estimador del error es menor que la tolerancia especificada. El proceso iterativo se considera completado.
!-------------------------------------------------------------------------------------------------------------------------------------
function isdbcg(n,b,itol,tol,iter,err,ierr,r,z,dz,rw,iw,bnrm,lenw,leniw)
use allocatacion
integer*4 ierr,iter,itol,n,lenw,leniw,iw(leniw),isdbcg
real*8 bnrm,err,tol,b(n),dz(n),r(n),rw(lenw),z(n),dnrm2

302 format (A7,I8,5X,A7,E11.4E2)

 !C�lculo del error estimado
 !--------------------------
 isdbcg=0
 if (itol.eq.1) then
   !err = ||residual||/||b|| (el s�mbolo ||^|| indica la norma 2 de ^).
   if (iter.eq.0) then
   bnrm=dnrm2(n,b)
   endif
   if (iter.gt.0) then
   err=dnrm2(n,r)/bnrm
   endif  
 elseif (itol.eq.2) then
   !err = ||inv(M)*residual||/||inv(M)*b||.
   if (iter.eq.0) then
     call dslui(n,b,dz,rw,iw,lenw,leniw)
     bnrm=dnrm2(n,dz)  
   endif
   err=dnrm2(n,z)/bnrm	
 else
   !Si se entra aqu�, itol no tiene un valor correcto.
   ierr=3
 endif

 !Escritura por pantalla de iteracion y del error.
 !------------------------------------------------
 if ((iter.eq.100).or.(iter.eq.500).or.(iter.eq.n)) then
 write(6,302) ' iter =',iter,'error =',err
 endif

 !Se eval�a si el error es mayor o menor que la tolerancia.
 !---------------------------------------------------------
 if (err.le.tol) then
 isdbcg=1
 endif
end

!--------------------------------------------------------------------------------------------------------------------------------------------
!Funci�n DNRM2
!Esta funci�n calcula la longitud Euclidiana (norma L2) de un vector (la raiz del cuadrado de sus componentes).
!'n' es el n�mero de elementos en el vector, 'dx' es el array donde est� el vector y tiene dimensi�n 'n' (los coeficientes no est�n 
!colocados espaciados). Es equivalente a la funci�n snrm.
!Se trabaja con dos contantes 'cutlo' y 'cuthi', y se eval�an cuatro fases.
!La fase 1 escanea componentes nulas. Se va a la fase 2 cuando una componente (distinta de cero) es menor o igual que 'cutlo'.
!Se va a la fase 3 cuando una componente es mayor que 'cutlo'. Se va a la fase 4 cuando una componente es mayor o igual que 'cuthi/n'.
!--------------------------------------------------------------------------------------------------------------------------------------------
function dnrm2 (n,dx)
integer*4 n,i,j  
real*8 dx(n),cutlo,cuthi,hitest,sum,xmax,dnrm2 
  
!Inicializaci�n previa de variables
!----------------------------------   
cutlo=sqrt(tiny(1.0_8)/epsilon(1.0_8))
cuthi=sqrt(huge(1.0_8))

if (n.le.0) then
  dnrm2=0.0_8 
  return
endif         			
sum=0.0_8 
	   
   !Comienzo del bucle principal
   !----------------------------
   i=1
   xmax=0.0_8
   !Fase 1. 'sum' es cero
   do j=i,n			   
    if (dx(j).eq.0.0_8) then	    
     if (j.eq.n) then   			   
	 dnrm2=xmax*sqrt(sum)
     !Sale de la subrutina con dnrm2=0 (no hay valores no nulos en el vector).
     return	
     endif
    else
    !Elimina los primeros ceros.	                            
    i=j
    exit
    endif
   enddo
	  
   !Eval�a para la primera componente no nula.
   if (abs(dx(i)).gt.cutlo) then	      
   !Establece hitest=cuthi/n.	  
   hitest=cuthi/n
   !Fase 3. 'sum' es de rango medio. Sin escala.
   do j=i,n	                         
    if (abs(dx(j)).ge.hitest) then
    !Se prepara para la fase 4.
    i=j
    sum=(sum/dx(i))/dx(i)
    xmax=abs(dx(i))
    sum=sum+(dx(i)/xmax)**2.0_8
    goto 200
    endif
    sum=sum+dx(j)**2.0_8
   enddo
   dnrm2=sqrt(sum)
   !Sale de la subrutina.
   return                             
   endif

   !Se prepara para la fase 2.
   !Primera componente no nula y menor que cutlo (valor muy peque�o).
   xmax=abs(dx(i))			         
   sum=sum+(dx(i)/xmax)**2.0_8
   i=i+1		                         
   if (i.le.n) then	                 
   !Fase 2. 'sum' es peque�o. Se escala para evitar underflow destructivo. 
70  if (abs(dx(i)).gt.cutlo) then     
    !Se prepara para la fase 3.
	sum=(sum*xmax)*xmax
    !Se establece hitest=cuthi/n                                         
    hitest=cuthi/n
    !Fase 3. 'sum' es de rango medio. Sin escala.
    do j=i,n                          
     if (abs(dx(j)).ge.hitest) then	   
     !Se prepara para la fase 4.
	 i=j
     sum=(sum/dx(i))/dx(i)
     xmax=abs(dx(i))
     sum=sum+(dx(i)/xmax)**2.0_8
     goto 200
     endif
    sum=sum+dx(j)**2.0_8
    enddo
    dnrm2=sqrt(sum)
    return
    endif
    !C�digo com�n para fases 2 y 4. En fase 4 'sum' es grande. Se escala para evitar overflow.
    if (abs(dx(i)).le.xmax) then
    sum=sum+(dx(i)/xmax)**2.0_8
    i=i+1	                         
     if (i.le.n) then                 
     goto 70
	 else
     !Fin del buble principal. C�lculo de la ra�z cuadrada y ajuste de la escala.		
 	 dnrm2=xmax*sqrt(sum)
	 return
	 endif
    endif
    sum=1.0_8+sum*(xmax/dx(i))**2.0_8 
    xmax=abs(dx(i))
    i=i+1                             
     if (i.le.n) then	             
	 goto 70 
     else
     !Fin del buble principal. C�lculo de la ra�z cuadrada y ajuste de la escala.
     dnrm2=xmax*sqrt(sum)
     return
     endif
   
   else
   !Fin del buble principal. C�lculo de la ra�z cuadrada y ajuste de la escala.
   dnrm2=xmax*sqrt(sum)
   return
   endif

200 do j=i+1,n                         
    !C�digo com�n para fases 2 y 4. En fase 4 sum es grande. Se escala para evitar overflow. 	   
	 if (abs(dx(j)).le.xmax) then	 
      sum=sum+(dx(j)/xmax)**2.0_8	 
      cycle                            
     endif
     sum=1.0_8+sum*(xmax/dx(j))**2.0_8	
     xmax=abs(dx(j))					                                
    enddo                                   
	     
!Fin del buble principal. C�lculo de la ra�z cuadrada y ajuste de la escala.
dnrm2=xmax*sqrt(sum)
end

!----------------------------------------------------------------------------------------------------- 
!FIN DE PROGRAMA
!-----------------------------------------------------------------------------------------------------
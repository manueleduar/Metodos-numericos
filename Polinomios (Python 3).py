'''
  Autor: Jose Villarreal

'''
# Clase Polinomio
class Polinomio(object):
  # Atributo de la clase polinomio
  #   Diccionario para guardar las constantes
  constantes = {}
  
  # Metodo de la clase polinomio
  #   Este metodo calcula el valor del polinomio en un punto
  #   x dado
  def evaluar(self,x):
    sum = 0
    print("x", "y")
    for k, v in self.constantes.items():
      sum += (int(x)** int(k) * int(v))
    print(x, sum)

# Funcion que crea y regresa un polinomio
def definirPolinomio():
  orden_poli = input("Orden del polinomio? ")
  print(orden_poli)
  if not orden_poli.isdigit():
    raise Exception("El orden debe de ser un entero")
  poli = Polinomio()
  for i in range(int(orden_poli)+1):
    const = input("Introduzca constante para x"+str(i) + " ")
    if not const.isdigit():
      raise Exception("La constante debe de ser un entero")
    poli.constantes[i] = const
  return poli
    
def main():
  p = None
  while True:
    if p == None:
      print("Presione la opcion deseada:")
      print("1 - Definir polinomio")
      print("0 - Salir")
      print(' ')
    else:
      print("Presione la opcion deseada:")
      print("1 - Definir polinomio")
      print("2 - Mostrar polinomio")
      print("3 - Evaluar polinomio")
      print("0 - Salir")
      print(' ')
    select = int(input())
    if select == 1:
      p = definirPolinomio()
    elif select == 0:
      break
    elif select == 2:
      for k, v in p.constantes.items():
        print("Constante del termino x^" + str(k) + " : "  + str(v))
    elif select == 3:
      x = ''
      while True:
        x = input('Que Valor de x? (s para terminar)')
        if x == 's':
          break
        p.evaluar(x)
main()
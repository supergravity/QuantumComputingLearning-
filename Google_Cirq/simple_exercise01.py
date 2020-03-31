# Simple Program in Cirq

# Import the Cirq Package 
import cirq

# Pick qubit
qubit = cirq.GridQubit(0,0)

#Create a circuit
circuit = cirq.Circuit.from_ops([
    cirq.X(qubit),  #NOT   
    cirq.measure(qubit, key='m') # Measurement
    ]
)

#Display the circuit
print("Circuit:")
print(circuit)

#Get a simulator to execute the circuit
simulator = cirq.Simulator()

#Simulate the circuit several times
result = simulator.run(circuit, repetitions = 10)

# Print the results
print("Results:")
print(result)


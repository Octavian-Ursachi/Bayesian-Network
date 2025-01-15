import tkinter as tk
import math
import xml.etree.ElementTree as ET
from tkinter import filedialog
from tkinter import messagebox

class Node:
    def __init__(self, name="default", x=0, y=0):
        self.x = x
        self.y = y
        self.name = name
        self.radius = 20
        self.inputs = []
        self.outputs = []
        self.prob_table = []
        self.observation = None

    def setRadius(self, radius):
        self.radius = radius

    def addInput(self, node):
        self.inputs.append(node)

    def addOutput(self, node):
        self.outputs.append(node)

class SimpleGraphApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Bayesian Network Viewer")

        # Frame for buttons
        self.button_frame = tk.Frame(root)
        self.button_frame.pack(pady=10)

        # Buttons
        self.import_tree_button = tk.Button(self.button_frame, text="Import Tree", command=self.import_tree)
        self.import_tree_button.grid(row=0, column=0, padx=5)

        self.modify_table_button = tk.Button(self.button_frame, text="Modify table", command=self.modify_table)
        self.modify_table_button.grid(row=0, column=1, padx=5)

        self.query_button = tk.Button(self.button_frame, text="Query", command=self.query)
        self.query_button.grid(row=0, column=2, padx=5)

        self.observation_button = tk.Button(self.button_frame, text="Make observation", command=self.observation)
        self.observation_button.grid(row=0, column=3, padx=5)

        self.action_buttons = [self.modify_table_button, self.query_button, self.observation_button]
        for button in self.action_buttons:
            button.config(bg="lightgray")
        self.selected_button = None

        # Canvas for drawing the graph
        self.canvas = tk.Canvas(root, width=600, height=600, bg="white")
        self.canvas.pack()

        # Initialize variables
        self.node_radius = 20
        self.treeNodes = []
        self.nodes_data = {}
        self.edges = []
        self.selected_node = None
        self.zoom_factor = 1.0  # Initial zoom level
        self.offset_x = 0  # X translation
        self.offset_y = 0  # Y translation

        # Bind mouse events
        self.canvas.bind("<Button-1>", self.on_click)
        self.canvas.bind("<B1-Motion>", self.on_drag)
        self.canvas.bind("<ButtonRelease-1>", self.on_release)
        self.canvas.bind("<MouseWheel>", self.on_zoom)  # Bind mouse wheel for zooming
        self.canvas.bind("<Double-1>", self.on_double_click)

    def initProbs(self, node):
        # Initialize all probabilities with 50%
        for i in range(0, pow(2, len(node.inputs)) * 2, 2):
            node.prob_table.append((0.5, 0.5))

    def update_button_state(self, active_button):
        # Reset all buttons to default style
        for button in self.action_buttons:
            if button == active_button:
                if self.selected_button == button:
                    self.selected_button = None
                    button.config(bg="lightgray")
                else:
                    self.selected_button = button
                    button.config(bg="lightblue")
            else:
                button.config(bg="lightgray")

    def import_tree(self):
        # Open file dialog to select XML file
        file_path = filedialog.askopenfilename(filetypes=[("XML Files", "*.xml")])
        self.treeNodes.clear()

        if file_path:
            self.parse_xml(file_path)
            self.normalize_positions()
            for node in self.treeNodes:
                self.initProbs(node)
            self.draw_graph(self.treeNodes)

    def modify_table(self):
        # Update the active button and its style
        self.update_button_state(self.modify_table_button)
        print("Modify table")

    def query(self):
        # Update the active button and its style
        self.update_button_state(self.query_button)
        print("Query")

    def observation(self):
        # Update the active button and its style
        self.update_button_state(self.observation_button)
        print("Make observation")

    def parse_xml(self, file_path):
        tree = ET.parse(file_path)
        root = tree.getroot()

        # Define namespaces to handle the XML namespaces correctly
        namespaces = {'': 'http://www.cs.ubc.ca/labs/lci/fopi/ve/XMLBIFv0_3'}

        # Extract VARIABLE nodes and their positions
        for variable in root.findall('.//VARIABLE', namespaces):
            name = variable.find('.//NAME', namespaces).text
            position = variable.find('.//PROPERTY', namespaces).text.split('=')[1].strip()[1:-1].split(',')
            x, y = float(position[0]), float(position[1])

            # Create node and add it to the list of tree nodes
            node = Node(name, x, y)
            self.treeNodes.append(node)
            self.nodes_data[name] = (x, y)

        # Extract the DEFINITION elements to determine dependencies (inputs and outputs)
        for definition in root.findall('.//DEFINITION', namespaces):
            for_node = definition.find('.//FOR', namespaces).text
            given_nodes = definition.findall('.//GIVEN', namespaces)
            for node in given_nodes:
                given_node_name = node.text
                # Find the nodes and link the relationships (edges)
                start_node = None
                end_node = None
                for node in self.treeNodes:
                    if node.name == for_node:
                        end_node = node
                    if node.name == given_node_name:
                        start_node = node
                if start_node and end_node:
                    start_node.addOutput(end_node)
                    end_node.addInput(start_node)

    def normalize_positions(self):
        # Get the min and max values for x and y
        min_x = min(node.x for node in self.treeNodes)
        max_x = max(node.x for node in self.treeNodes)
        min_y = min(node.y for node in self.treeNodes)
        max_y = max(node.y for node in self.treeNodes)

        # Normalize the positions to fit within the canvas size
        for node in self.treeNodes:
            normalized_x = (node.x - min_x) / (max_x - min_x) * 600  # canvas width
            normalized_y = (node.y - min_y) / (max_y - min_y) * 600  # canvas height
            node.x = normalized_x
            node.y = normalized_y

    def draw_graph(self, nodes):
        # Clear the canvas before drawing new graph
        self.canvas.delete("all")

        # Draw edges
        for node in nodes:
            for output_node in node.outputs:
                # Get coordinates of start and end nodes
                x1, y1 = node.x, node.y
                x2, y2 = output_node.x, output_node.y

                # Adjust line length to avoid overlapping the node
                x1, y1, x2, y2 = self.adjust_line(x1, y1, x2, y2, self.node_radius)

                # Draw the line (edge)
                self.canvas.create_line(
                    x1 * self.zoom_factor + self.offset_x, y1 * self.zoom_factor + self.offset_y,
                    x2 * self.zoom_factor + self.offset_x, y2 * self.zoom_factor + self.offset_y,
                    width=2
                )

                # Draw a triangle (arrowhead) at the end of the line
                self.draw_arrowhead(x1 * self.zoom_factor + self.offset_x, y1 * self.zoom_factor + self.offset_y,
                                    x2 * self.zoom_factor + self.offset_x, y2 * self.zoom_factor + self.offset_y)

        # Draw nodes
        for node in nodes:
            x = node.x
            y = node.y
            # Draw node (circle)
            self.canvas.create_oval(
                (x - self.node_radius) * self.zoom_factor + self.offset_x,
                (y - self.node_radius) * self.zoom_factor + self.offset_y,
                (x + self.node_radius) * self.zoom_factor + self.offset_x,
                (y + self.node_radius) * self.zoom_factor + self.offset_y,
                fill="lightblue", tags=node.name
            )
            # Draw node label
            self.canvas.create_text(
                x * self.zoom_factor + self.offset_x,
                y * self.zoom_factor + self.offset_y,
                text=node.name, font=("Arial", 12, "bold")
            )

    def adjust_line(self, x1, y1, x2, y2, offset):
        """Adjust the line to stop at the edge of the node."""
        angle = math.atan2(y2 - y1, x2 - x1)
        x2 -= offset * math.cos(angle)
        y2 -= offset * math.sin(angle)
        return x1, y1, x2, y2

    def draw_arrowhead(self, x1, y1, x2, y2):
        # Calculate the direction vector for the line
        angle = math.atan2(y2 - y1, x2 - x1)

        # Define arrowhead size
        arrow_length = 10
        arrow_width = 5

        # Calculate the base of the triangle
        tip_x, tip_y = x2, y2
        base_left_x = x2 - arrow_length * math.cos(angle) + arrow_width * math.sin(angle)
        base_left_y = y2 - arrow_length * math.sin(angle) - arrow_width * math.cos(angle)
        base_right_x = x2 - arrow_length * math.cos(angle) - arrow_width * math.sin(angle)
        base_right_y = y2 - arrow_length * math.sin(angle) + arrow_width * math.cos(angle)

        # Draw the triangle
        self.canvas.create_polygon(
            tip_x, tip_y,
            base_left_x, base_left_y,
            base_right_x, base_right_y,
            fill="black"
        )

    def open_prob_table_window(self, node):
        # Create a new top-level window
        new_window = tk.Toplevel(self.root)
        new_window.title("New Window")

        # Block interaction with the main window
        new_window.grab_set()

        # Table frame 
        nr_inputs = len(node.inputs)

        # 2^x input nodes
        heigh = pow(2, nr_inputs)
        # input nodes + 2 entry
        width = nr_inputs + 2

        table_frame = tk.Frame(new_window)
        table_frame.grid(row = 0, column=0,pady=10)

        labels = []
        for nod in node.inputs:
            label = tk.Label(table_frame, text=f'{nod.name}')
            labels.append(label)
        labels.append(tk.Label(table_frame, text=f'P ( {node.name} = T )'))
        labels.append(tk.Label(table_frame, text=f'P ( {node.name} = F )'))

        entries = []
        for i in range(heigh * 2):
            entries.append(tk.Entry(table_frame, text=""))
             # Pre-fill the entry with data from node.prob_table if available
            if node.prob_table and len(node.prob_table) > (i // 2):
                entry_value = node.prob_table[i // 2][i % 2]
                entries[i].insert(0, str(entry_value))

    
        tf = 0
        nr_bits = f'0{nr_inputs}b'
        for i in range(heigh + 1):
            for j in range(width):
                if i == 0:
                    labels[j].grid(row = i, column = j)
                else:
                    if j > nr_inputs - 1:
                        entry_index = (i - 1) * 2 + (j - nr_inputs)
                        entries[entry_index].grid(row = i, column = j)
                    else:
                        if i > 0 :
                            binary_rep = format(tf, nr_bits)
                            if binary_rep[j] == '0':
                                tk.Label(table_frame, text='T').grid(row = i, column= j)
                            else:
                                tk.Label(table_frame, text='F').grid(row = i, column= j)
            if i > 0:
                tf += 1

        # Create a frame for buttons
        button_frame = tk.Frame(new_window)
        button_frame.grid(row=1, column=0, pady=10)  # Grid it, not pack it

        # Buttons in the frame
        ok_button = tk.Button(button_frame, text="OK", command=lambda: self.table_ok(entries, node, new_window))
        ok_button.grid(row=0, column=0, padx=5)

        close_button = tk.Button(button_frame, text="Close", command=new_window.destroy)
        close_button.grid(row=0, column=1, padx=5)

    def table_ok(self, entries, node, window):
        node.prob_table.clear()
        for i in range(0, len(entries), 2):
            value1 = entries[i].get().strip()  # strip to avoid leading/trailing spaces
            value2 = entries[i + 1].get().strip()
            
            try:
                value1 = float(value1)
                value2 = float(value2)
            except ValueError:
                messagebox.showerror("Invalid Input", f"Invalid input: {value1}, {value2}. Please enter numeric values.")
                return  # Stop further processing

            if value1 + value2 != 1:
                messagebox.showwarning("Invalid Sum", f"The sum of {value1} and {value2} is not equal to 1. Please correct it.")
                return  # Stop further processing
            
            node.prob_table.append((value1, value2))
        
        window.destroy()

    def compute_node_probability(self, node):
        """
        Calculate the probability of a node using the Bayesian network logic.
        
        :param node: The Node object for which to calculate probability.
        :return: A dictionary with the probabilities {True: P(Node=True), False: P(Node=False)}.
        """
        # If the node has no parents, use marginal probabilities
        if not node.inputs:
            return {True: node.prob_table[0][0], False: node.prob_table[0][1]}

        # Initialize probabilities
        probabilities = {True: 0.0, False: 0.0}

        # Generate all possible parent states
        parent_states = self.generate_parent_states(node.inputs)

        # Compute probabilities for each state
        for idx, state in enumerate(parent_states):
            # Recursive probability calculation for parent states
            prob_state = self.calculate_state_probability(state, node.inputs)
            probabilities[True] += prob_state * node.prob_table[idx][0]
            probabilities[False] += prob_state * node.prob_table[idx][1]

        return probabilities

    def generate_parent_states(self, parents):
        """
        Generate all possible states for parent nodes.
        
        :param parents: List of parent nodes.
        :return: List of dictionaries with all combinations of parent states.
        """
        from itertools import product
        parent_names = [parent.name for parent in parents]
        parent_values = [[True, False] for _ in parents]
        combinations = product(*parent_values)
        return [dict(zip(parent_names, values)) for values in combinations]

    def calculate_state_probability(self, state, parents):
        """
        Calculate the joint probability of a state for the parent nodes.
        
        :param state: A dictionary representing a combination of parent states.
        :param parents: List of parent nodes.
        :return: Joint probability of the given state.
        """
        probability = 1.0
        for parent in parents:
            # Calculate probability of the parent recursively
            parent_prob = self.compute_node_probability(parent)
            probability *= parent_prob[True] if state[parent.name] else parent_prob[False]
        return probability

    def open_query_popup(self, node):
        prob = self.compute_node_probability(node)

        # Create a new popup window
        popup = tk.Toplevel(self.root)
        popup.title(f"Probabilities for {node.name}")

        # Display probabilities in the popup
        frame = tk.Frame(popup)
        frame.pack(padx=10, pady=10)

        tk.Label(frame, text=f"Probability of {node.name} = True:").grid(row=0, column=0, sticky="w")
        tk.Label(frame, text=f"{prob[True]:.4f}").grid(row=0, column=1, sticky="e")

        tk.Label(frame, text=f"Probability of {node.name} = False:").grid(row=1, column=0, sticky="w")
        tk.Label(frame, text=f"{prob[False]:.4f}").grid(row=1, column=1, sticky="e")

        # Add a close button
        close_button = tk.Button(frame, text="Close", command=popup.destroy)
        close_button.grid(row=2, column=0, columnspan=2, pady=10)

    def open_observation_window(self, node):
        new_window = tk.Toplevel(self.root)
        new_window.title("Observation Window")

        new_window.grab_set()

        frame = tk.Frame(new_window)
        frame.grid(row=0, column=0, padx=10, pady=10)

        label = tk.Label(frame, text="Select an option:")
        label.grid(row=0, column=0, pady=5)

        options = ["None", "True", "False"]
        selected_option = tk.StringVar(value=options[0])

        option_menu = tk.OptionMenu(frame, selected_option, *options)
        option_menu.grid(row=1, column=0, padx=5, pady=5)

        ok_button = tk.Button(frame, text="OK", command=lambda: self.observation_ok(node, new_window, selected_option.get()))
        ok_button.grid(row=2, column=0, pady=10)

    def observation_ok(self, node, window, selected_value):
        node.observation = selected_value
        window.destroy()

    def open_new_window(self, node):
        if self.selected_button == self.modify_table_button:
            self.open_prob_table_window(node)
        elif self.selected_button == self.query_button:
            self.open_query_popup(node)
        elif self.selected_button == self.observation_button:
            self.open_observation_window(node)

    def on_double_click(self, event):
        # Normalize the coordinates with the current zoom and offset
        click_x = (event.x - self.offset_x) / self.zoom_factor
        click_y = (event.y - self.offset_y) / self.zoom_factor

        # Check if any node is clicked
        for node in self.treeNodes:
            x1, y1 = node.x - self.node_radius, node.y - self.node_radius
            x2, y2 = node.x + self.node_radius, node.y + self.node_radius
            if x1 <= click_x <= x2 and y1 <= click_y <= y2:
                self.open_new_window(node)
                break

    def on_click(self, event):
        # Scale the mouse position to account for zoom
        click_x = (event.x - self.offset_x) / self.zoom_factor
        click_y = (event.y - self.offset_y) / self.zoom_factor

        # Check if any node is clicked
        for node in self.treeNodes:
            x1, y1 = node.x - self.node_radius, node.y - self.node_radius
            x2, y2 = node.x + self.node_radius, node.y + self.node_radius
            if x1 <= click_x <= x2 and y1 <= click_y <= y2:
                self.selected_node = node
                break

    def on_drag(self, event):
        # If a node is selected, move it with the mouse
        if self.selected_node:
            # Scale the drag position to account for zoom
            drag_x = (event.x - self.offset_x) / self.zoom_factor
            drag_y = (event.y - self.offset_y) / self.zoom_factor
            self.selected_node.x = drag_x
            self.selected_node.y = drag_y
            self.draw_graph(self.treeNodes)

    def on_release(self, event):
        # Release the node when mouse button is released
        self.selected_node = None

    def on_zoom(self, event):
        """Handle zoom in/out using the mouse wheel"""
        # Zoom in (up) or zoom out (down)
        zoom_step = 0.1
        old_zoom_factor = self.zoom_factor

        if event.delta > 0:
            self.zoom_factor += zoom_step  # Zoom in
        else:
            self.zoom_factor -= zoom_step  # Zoom out

        self.zoom_factor = max(0.1, self.zoom_factor)  # Prevent zooming out too much

        # Scale the positions of all nodes based on the new zoom factor
        scale_factor = self.zoom_factor / old_zoom_factor

        for node in self.treeNodes:
            node.x *= scale_factor
            node.y *= scale_factor

        self.offset_x = event.x - scale_factor * (event.x - self.offset_x)
        self.offset_y = event.y - scale_factor * (event.y - self.offset_y)

        self.draw_graph(self.treeNodes)  # Redraw the graph with the updated zoom level

if __name__ == "__main__":
    root = tk.Tk()
    app = SimpleGraphApp(root)
    root.mainloop()
